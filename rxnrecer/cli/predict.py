import sys
import os
import time
from types import SimpleNamespace
import argparse
from collections import defaultdict
import tempfile
from typing import Any
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import torch

from rxnrecer.utils import file_utils as ftool
from rxnrecer.utils import bio_utils as butils
from rxnrecer.lib.llm import qa as llm_qa
from rxnrecer.lib.model import mactive as Mactive
from rxnrecer.utils import format_utils

if __package__ in {None, ""}:
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

from rxnrecer.config import config as cfg

def load_model(model_weight_path=None, device=None):
    
    
    from rxnrecer.lib.model import mactive as Mactive

    if model_weight_path is None:
        model_weight_path = cfg.FILE_MOLEL_PRODUCTION_BEST_MODEL

    if device is None:
        resolved_device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        resolved_device = torch.device(device)

    mcfg = SimpleNamespace(
    # Model parameters
    batch_size=1,
    esm_out_dim=1280,
    gru_h_dim=512,
    att_dim=32,
    dropout_rate=0.2,
    freeze_esm_layers = 32, # Number of frozen layers
    output_dimensions=10479,
    device = resolved_device,
    model_weight_path=model_weight_path,
    dict_path = cfg.FILE_DS_DICT_ID2RXN
    )

    print(f'Use device: {mcfg.device}')
    model = Mactive.BGRU(
        input_dimensions=mcfg.esm_out_dim,
        device = mcfg.device,
        gru_h_size=mcfg.gru_h_dim,
        attention_size=mcfg.att_dim,
        dropout=mcfg.dropout_rate,
        output_dimensions=mcfg.output_dimensions,
        freeze_esm_layers=mcfg.freeze_esm_layers,
    )
    
    return model, mcfg



# region 1. Load computation data
def load_data(input_data) -> pd.DataFrame:


    """
    Load computation data, supports two input formats:
      1) Pass FASTA file path (str) -> auto parse to DataFrame
      2) Pass DataFrame directly -> must contain 'seq' and 'uniprot_id' columns
    """
    if isinstance(input_data, str):
        # Assume input is FASTA file path
        print(f'Detected input is a FASTA file: {input_data}')
        input_df = ftool.fasta_to_dataframe(fasta_file=input_data)
    elif isinstance(input_data, pd.DataFrame):
        input_df = input_data.copy().reset_index(drop=True)
    else:
        raise ValueError("input_data should be either a path to FASTA or a pandas DataFrame.")
    
    # Ensure required columns exist
    if 'seq' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'seq'.")
    if 'uniprot_id' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'uniprot_id'.")
    
    return input_df
# endregion 




#region 2. Save results
def normalize_output_file(output_file, output_format) -> str:
    """Return an output path whose suffix matches the selected format."""
    expected_suffix = f".{output_format}"
    root, suffix = os.path.splitext(output_file)

    if suffix.lower() == expected_suffix:
        return output_file

    if suffix.lower() in {".tsv", ".json"}:
        return root + expected_suffix

    return output_file + expected_suffix


def save_data(resdf: pd.DataFrame, output_file: str, output_format='tsv') -> None:

    """
    Save results to file (TSV/JSON)
    resdf: Result DataFrame
    output_file: Output file path
    output_format: 'tsv' or 'json'
    """
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    def _atomic_write(writer, suffix):
        fd, temp_path = tempfile.mkstemp(dir=output_dir or None, suffix=suffix)
        os.close(fd)
        try:
            writer(temp_path)
            os.replace(temp_path, output_file)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def format_cell(x):
        return butils.format_obj(x, 4)

    resdf = resdf.applymap(format_cell)
    if output_format == 'tsv':
        resdf[["reaction_ec", "reaction_equation", "reaction_equation_ref_chebi"]] = resdf['rxn_details'].apply(butils.simplify_rxn_details_fields).apply(pd.Series)
        resdf = resdf[['input_id', 'RXNRECer', 'RXNRECer_with_prob', 'reaction_ec', 'reaction_equation', 'reaction_equation_ref_chebi']]
        def _write_tsv(path):
            with open(path, "w", encoding="utf-8") as handle:
                resdf.to_csv(handle, index=False, sep='\t', float_format="%.4f")
                handle.flush()
                os.fsync(handle.fileno())
        _atomic_write(_write_tsv, ".tsv")
    elif output_format == 'json':
        def _write_json(path):
            with open(path, "w", encoding="utf-8") as handle:
                resdf.to_json(handle, orient='records', indent=4)
                handle.flush()
                os.fsync(handle.fileno())
        _atomic_write(_write_json, ".json")
    else:
        print(f'Error: Invalid output format {output_format}. Skip saving.')
# endregion 


def _split_chunks_even(indexed_df: pd.DataFrame, num_chunks: int) -> list[pd.DataFrame]:
    return [chunk.copy() for chunk in np.array_split(indexed_df, num_chunks) if not chunk.empty]


def _split_chunks_by_sequence_length(indexed_df: pd.DataFrame, num_chunks: int) -> list[pd.DataFrame]:
    work_df = indexed_df.copy()
    work_df["_seq_len"] = work_df["seq"].fillna("").astype(str).str.len()
    work_df = work_df.sort_values("_seq_len", ascending=False).reset_index(drop=True)
    buckets = [[] for _ in range(num_chunks)]
    bucket_lengths = [0] * num_chunks

    for _, row in work_df.iterrows():
        bucket_idx = min(range(num_chunks), key=lambda idx: bucket_lengths[idx])
        buckets[bucket_idx].append(row.to_dict())
        bucket_lengths[bucket_idx] += int(row["_seq_len"])

    chunks = []
    for bucket in buckets:
        if not bucket:
            continue
        chunk_df = pd.DataFrame(bucket).drop(columns=["_seq_len"])
        chunks.append(chunk_df.reset_index(drop=True))
    return chunks


def _build_prediction_chunks(input_df: pd.DataFrame, num_chunks: int, split_mode: str = "balanced_length") -> list[pd.DataFrame]:
    indexed_df = input_df.reset_index(drop=False).rename(columns={"index": "_orig_idx"})
    if split_mode == "even":
        return _split_chunks_even(indexed_df, num_chunks)
    return _split_chunks_by_sequence_length(indexed_df, num_chunks)


def _print_chunk_stats(chunks: list[pd.DataFrame], devices: list[str], split_mode: str) -> None:
    print(f"[predict] split_mode={split_mode}")
    for device_name, chunk in zip(devices, chunks):
        seq_len = chunk["seq"].fillna("").astype(str).str.len()
        print(
            f"[predict] {device_name}: samples={len(chunk)} "
            f"total_len={int(seq_len.sum())} max_len={int(seq_len.max()) if len(seq_len) else 0}"
        )


def _predict_chunk_worker(
    input_df: pd.DataFrame,
    model_weight_path: str,
    mode: str,
    batch_size: int,
    rxn_details: bool,
    device: str,
    progress_desc: str | None = None,
    progress_position: int | None = None,
) -> pd.DataFrame:
    model, mcfg = load_model(model_weight_path=model_weight_path, device=device)
    return single_batch_run_prediction(
        input_df=input_df,
        model=model,
        mcfg=mcfg,
        mode=mode,
        batch_size=batch_size,
        rxn_details=rxn_details,
        show_progress=True,
        progress_desc=progress_desc,
        progress_position=progress_position,
        devices=None,
    ).drop_duplicates(subset=["input_id"], keep="first")


def step_by_step_prediction(input_data, mode='s1', batch_size=100,  rxn_details=True, devices=None, split_mode="balanced_length") -> pd.DataFrame:
    if mode == 's3':
        if not cfg.LLM_API_KEY or not cfg.LLM_API_URL:
            print("Error: LLM API key and URL are required for S3 mode!")
            return pd.DataFrame()
    
    input_df = load_data(input_data)
    print(f'Step 1: Preparing input data, loading {len(input_df)} proteins')

    # Step 2: Load predictive model
    print('Step 2: Loading predictive model')
    model, mcfg = load_model(model_weight_path=cfg.FILE_MOLEL_PRODUCTION_BEST_MODEL) 

    print(f'Step 3: Running prediction on {len(input_df)} proteins')
    if input_df.empty:
        return pd.DataFrame()

    fres = single_batch_run_prediction(
        input_df=input_df,
        model=model,
        mcfg=mcfg,
        mode=mode,
        batch_size=batch_size,
        rxn_details=rxn_details,
        devices=devices,
        split_mode=split_mode,
    )
    fres = fres.drop_duplicates(subset=['input_id'], keep='first')
    return fres
    


#region RXNRECer Prediction API
def single_batch_run_prediction(
    input_df: pd.DataFrame,
    model: Any,
    mcfg: SimpleNamespace,
    mode='s1',
    batch_size=2,
    rxn_details=True,
    show_progress=True,
    progress_desc=None,
    progress_position=None,
    devices=None,
    split_mode="balanced_length",
) -> pd.DataFrame:

    """
    Perform RXNRECer prediction on input DataFrame and return prediction results.
    Args:
      - input_df: Input DataFrame with 'seq' and 'uniprot_id' columns
      - model: Pre-trained model
      - mcfg: Model configuration
      - mode: Prediction mode ('s1' or 's2')
    """
    input_df = input_df.reset_index(drop=True)
    if devices is not None and not isinstance(devices, str):
        device_list = [item for item in devices if item]
        if len(device_list) > 1:
            chunks = _build_prediction_chunks(input_df=input_df, num_chunks=len(device_list), split_mode=split_mode)
            _print_chunk_stats(chunks, device_list, split_mode)

            with ProcessPoolExecutor(max_workers=len(chunks), mp_context=get_context("spawn")) as executor:
                futures = [
                    executor.submit(
                        _predict_chunk_worker,
                        chunk[["uniprot_id", "seq"]],
                        mcfg.model_weight_path,
                        mode,
                        batch_size,
                        rxn_details,
                        device_name,
                        f"Predicting on {device_name}",
                        worker_idx,
                    )
                    for worker_idx, (chunk, device_name) in enumerate(zip(chunks, device_list))
                ]

                result_frames = []
                for chunk, future in zip(chunks, futures):
                    result_frame = future.result()
                    result_frame = result_frame.merge(
                        chunk[["_orig_idx", "uniprot_id"]].rename(columns={"uniprot_id": "input_id"}),
                        on="input_id",
                        how="left",
                    )
                    result_frames.append(result_frame)

            return (
                pd.concat(result_frames, axis=0, ignore_index=True)
                .sort_values("_orig_idx")
                .drop(columns=["_orig_idx"])
                .reset_index(drop=True)
            )
    # RXNRECer-S1
    # print('Running RXNRECer-S1 ...')
    res, res_prob = Mactive.predict_sequences(
        model=model,
        sequences=input_df.seq,
        model_weight_path=mcfg.model_weight_path,
        dict_path=mcfg.dict_path,
        batch_size=batch_size,
        device=mcfg.device,
        show_progress=show_progress,
        progress_desc=progress_desc,
        progress_position=progress_position,
    )
    # Integrate prediction results
    res_df_s1 = input_df[['uniprot_id']].reset_index(drop=True).copy()
    res_df_s1['RXNRECer'] = res
    res_df_s1['RXNRECer_with_prob'] = res_prob
    res_df_s1 = res_df_s1.rename(columns={'uniprot_id': 'input_id'})
    
    res_df_s1 = refine_prediction_table(res_df_s1)
    
    if mode == 's1':
        # Get reaction details
        rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
        res_df_s1= butils.get_rxn_details_batch(df_rxns=res_df_s1, rxn_bank=rxn_bank, rxn_id_column='RXNRECer')
        if rxn_details:
            res_df_s1['rxn_details'] = res_df_s1.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, RXN_details=x.RXN_details, mode='s2'), axis=1).tolist()
            res_df_s1 = res_df_s1[['input_id','RXNRECer', 'RXNRECer_with_prob', 'rxn_details']]
        else:
            res_df_s1 = res_df_s1[['input_id','RXNRECer', 'RXNRECer_with_prob']]
        return res_df_s1


    # Ensemble learning
    if mode == 's2':
        print('Running RXNRECer-S2 ...')
        res_df_s2 = get_ensemble(input_df=input_df[['uniprot_id', 'seq']], rxnrecer_df=res_df_s1)  # Ensemble learning
        # Get reaction details
        rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
        res_df_s2= butils.get_rxn_details_batch(df_rxns=res_df_s2, rxn_bank=rxn_bank, rxn_id_column='RXNRECer')
        # Format output
        if rxn_details:
            res_df_s2['rxn_details'] = res_df_s2.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, RXN_details=x.RXN_details, mode='s2'), axis=1).tolist()
            res_df_s2 = res_df_s2[['input_id','RXNRECer', 'RXNRECer_with_prob', 'rxn_details']]
        else:
            res_df_s2 = res_df_s2[['input_id','RXNRECer', 'RXNRECer_with_prob']]
        return res_df_s2
    
    if mode == 's3':
        print('Running RXNRECer-S3 ...')
        res_df_s2 = get_ensemble(input_df=input_df[['uniprot_id', 'seq']], rxnrecer_df=res_df_s1).rename(columns={'RXNRECer': 'RXNRECer-S2'})
        s3_input_df = res_df_s2.merge(res_df_s1[['input_id', 'RXNRECer']].rename(columns={'RXNRECer': 'RXNRECer-S1'}), on='input_id', how='left'
                            ).merge(input_df.rename(columns={'uniprot_id': 'input_id'}), on='input_id', how='left')
        res_df_s3 = llm_qa.batch_chat(res_rxnrecer=s3_input_df, api_key=cfg.LLM_API_KEY, api_url=cfg.LLM_API_URL, llm_model=cfg.LLM_MODEL)
        res_df_s3['rxn_details'] = res_df_s3.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, 
                                                                                            RXNRECER_S3=x['RXNRECER-S3'], 
                                                                                            RXN_details=x.RXN_details, mode='s3'), axis=1)
        res_df_s3 = res_df_s3[['input_id','RXNRECer-S2', 'RXNRECer_with_prob', 'rxn_details']].rename(columns={'RXNRECer-S2': 'RXNRECer'})
        
        return res_df_s3
        
# endregion

#region Handle enzyme/non-enzyme mixed cases in ensemble
def res_refinement(rxn_prob):
    """
    Refine enzyme prediction results with rules:
    1. Single result (enzyme or non-enzyme) -> keep as is
    2. If '-' has highest probability -> non-enzyme, keep '-'
    3. If '-' has lower probability -> remove '-'
    """
    from rxnrecer.config import config as cfg

    if len(rxn_prob) == 1:
        return cfg.SPLITER.join(rxn_prob.keys()), rxn_prob

    sorted_items = sorted(rxn_prob.items(), key=lambda x: x[1], reverse=True)

    # If '-' has highest probability -> non-enzyme, keep
    if sorted_items[0][0] == '-':
        return '-', {'-': rxn_prob['-']}

    # If '-' exists and not highest -> remove regardless of position
    if '-' in rxn_prob:
        rxn_prob.pop('-')

    return cfg.SPLITER.join(rxn_prob.keys()), rxn_prob


def refine_prediction_table(
    res_df: pd.DataFrame,
    rxn_col: str = "RXNRECer",
    prob_col: str = "RXNRECer_with_prob",
) -> pd.DataFrame:
    """统一清理预测结果里的酶/非酶冲突。

    允许的合法状态只有两种：
    1. 只有 `-`
    2. 一个或多个真实反应，不包含 `-`

    因此如果同一条预测里同时出现 `-` 和真实反应：
    - 若 `-` 概率最高，则整条视为非酶，只保留 `-`
    - 否则删除 `-`，保留真实反应
    """
    working = res_df.copy()
    working[[rxn_col, prob_col]] = working.apply(
        lambda row: pd.Series(res_refinement(dict(row[prob_col]))),
        axis=1,
    )
    return working
# endregion


def get_ensemble(input_df, rxnrecer_df):
    import pandas as pd
    from rxnrecer.config import config as cfg
    from rxnrecer.lib.ml import mlpredict as predRXN

    """
    Ensemble multiple methods for protein reaction prediction and handle enzyme/non-enzyme mixed cases.

    Args:
    - input_df: Protein info DataFrame with 'uniprot_id' column
    - rxnrecer_df: RXNRECer prediction results with 'input_id' and 'RXNRECer_with_prob' columns

    Returns:
    - res_df: DataFrame with input_id, ensemble prediction string (RXNRECer) and probability dict (RXNRECer_with_prob)
    """

    def _fallback_rxn_frame(column_name):
        return input_df[["uniprot_id"]].assign(**{column_name: "NO-PREDICTION"})

    def _fallback_similarity_frame():
        empty_prediction = [[("-", 0.0)]]
        frame = input_df[["uniprot_id"]].copy()
        frame["esm"] = empty_prediction * len(frame)
        frame["t5"] = empty_prediction * len(frame)
        return frame

    def _run_or_fallback(label, fn, fallback_factory):
        try:
            return fn()
        except Exception as exc:
            print(f"Warning: {label} failed, fallback to NO-PREDICTION. Reason: {exc}")
            return fallback_factory()

    # Call various prediction methods (MSA, CatFam, ECRECer, T5, RXNRECer)
    res_msa = _run_or_fallback("MSA", lambda: predRXN.getmsa(df_test=input_df, k=1), lambda: _fallback_rxn_frame("rxn_msa"))
    res_catfam = _run_or_fallback("CatFam", lambda: predRXN.getcatfam(df_test=input_df), lambda: _fallback_rxn_frame("rxn_catfam"))
    res_ecrecer = _run_or_fallback("ECRECer", lambda: predRXN.getecrecer(df_test=input_df), lambda: _fallback_rxn_frame("rxn_ecrecer"))
    res_t5 = _run_or_fallback("T5", lambda: predRXN.getT5(df_test=input_df, topk=1), _fallback_similarity_frame)
    res_rxnrecer = rxnrecer_df.copy().rename(columns={'input_id': 'uniprot_id'})

    # Merge different model results into same DataFrame (left join on uniprot_id)
    baggingdf = res_rxnrecer.merge(
                    res_msa[['uniprot_id', 'rxn_msa']], on='uniprot_id', how='left'
                ).merge(
                    res_catfam[['uniprot_id', 'rxn_catfam']], on='uniprot_id', how='left'
                ).merge(
                    res_ecrecer[['uniprot_id', 'rxn_ecrecer']], on='uniprot_id', how='left'
                ).merge(
                    res_t5, on='uniprot_id', how='left'
                )

    # Standardize missing predictions: convert NO-PREDICTION, None and NaN to '-'
    baggingdf.replace(['NO-PREDICTION', 'None'], '-', inplace=True)
    baggingdf.fillna('-', inplace=True)

    # Call ensemble method for fusion, return reaction prediction dict with probabilities for each protein
    baggingdf['ensemble'] = baggingdf.apply(lambda x: integrateEnsemble(
        esm=x.esm,
        t5=x.t5,
        rxn_recer=x.RXNRECer_with_prob,
        rxn_msa=x.rxn_msa,
        rxn_catfam=x.rxn_catfam,
        rxn_ecrecer=x.rxn_ecrecer
    ), axis=1)

    # Extract ensemble results into two columns:
    # - RXNRECer: concatenated reaction ids (string)
    # - RXNRECer_with_prob: original ensemble dict
    baggingdf['RXNRECer_with_prob'] = baggingdf['ensemble']
    baggingdf['RXNRECer'] = baggingdf['ensemble'].apply(lambda x: cfg.SPLITER.join(x.keys()))

    # Restore input_id naming for consistency with external interface
    baggingdf.rename(columns={'uniprot_id': 'input_id'}, inplace=True)

    # Keep only required output fields
    res_df = baggingdf[['input_id', 'RXNRECer', 'RXNRECer_with_prob']].copy()

    return refine_prediction_table(res_df)
    

def integrateEnsemble(esm, t5, rxn_recer, rxn_msa, rxn_catfam, rxn_ecrecer):
    """
    整合不同来源的 RHEA ID 概率，优先保留较高的概率值，并对单个 RHEA ID 进行归并。
    
    参数:
        esm: 包含 (RHEA_IDs, 概率) 的元组列表，这里取第一个元素
        t5: 包含 (RHEA_IDs, 概率) 的元组列表，这里取第一个元素
        rxn_recer: 字典，键为 RHEA_ID，值为对应的概率
        rxn_msa: 字符串，MSA 方法产生的 RHEA ID（可能包含多个，用分号分隔）
        rxn_catfam: 字符串，分类家族方法产生的 RHEA ID（可能包含多个，用分号分隔）
        rxn_ecrecer: 字符串，EC 方法产生的 RHEA ID（可能包含多个，用分号分隔）
    
    返回:
        dict: 以单个 RHEA ID 为键，对应最大概率为值的字典，按概率降序排序
    """
    # -------------------------------
    # 1. 聚合 esm 与 t5 的预测结果
    # -------------------------------
    # 使用 defaultdict(list) 将同一组（可能包含多个 RHEA ID，以分号分隔）的概率收集在一起
    aggregated = defaultdict(list)
    for prediction in (esm[0], t5[0]):
        # 如果 RHEA ID 为 None，则替换为 '-' 以保持一致性
        rhea_ids = prediction[0] if prediction[0] is not None else '-'
        aggregated[rhea_ids].append(prediction[1])
    
    # -------------------------------
    # 2. 添加 rxn_recer 的预测结果
    # -------------------------------
    # 同样确保 None 替换为 '-'
    for rhea_id, prob in rxn_recer.items():
        key = rhea_id if rhea_id is not None else '-'
        aggregated[key].append(prob)
    
    # -------------------------------
    # 3. 计算每组 RHEA ID 的最高概率
    # -------------------------------
    # aggregated 的键可能包含多个 RHEA ID（如 "ID1;ID2"），此处对每组取最大概率
    group_max_prob = {group: max(probs) for group, probs in aggregated.items()}
    
    # -------------------------------
    # 4. 将 RHEA ID 组拆分成单个 RHEA ID，并保留较高概率
    # -------------------------------
    direct_prob = {}
    for group, prob in group_max_prob.items():
        # 拆分可能由分号连接的多个 RHEA ID
        for rhea_id in group.split(';'):
            # 如果同一 RHEA ID 来自多个组，保留概率较大的那个
            direct_prob[rhea_id] = max(direct_prob.get(rhea_id, 0), prob)
    
    # -------------------------------
    # 5. 构建 EC 方法手动添加的 RHEA ID 字典
    # -------------------------------
    # 将 rxn_msa, rxn_catfam, rxn_ecrecer 三个字符串拼接后，以分号拆分，得到唯一的 RHEA ID 集合
    ec_ids = set(f"{rxn_msa};{rxn_catfam};{rxn_ecrecer}".split(';'))
    # 对每个 EC 方法产生的 RHEA ID 赋予固定概率 0.7777
    ec_prob = {rhea_id: 0.7777 for rhea_id in ec_ids}
    
    # -------------------------------
    # 6. 合并来自 direct_prob 与 ec_prob 的结果
    # -------------------------------
    # 对于相同的 RHEA ID，取两边概率中的较大值
    merged_probs = {}
    for mapping in (ec_prob, direct_prob):
        for rhea_id, prob in mapping.items():
            merged_probs[rhea_id] = max(merged_probs.get(rhea_id, 0), prob)
    
    # -------------------------------
    # 7. 按概率降序排序并返回结果
    # -------------------------------
    sorted_result = dict(sorted(merged_probs.items(), key=lambda item: item[1], reverse=True))
    return sorted_result

    

def main():
    """Main function for command line interface"""
    from rxnrecer.config import config as cfg
    from rxnrecer.utils import file_utils as ftool

    start_time = time.perf_counter()
    
    # 1. 解析命令行参数
    parser = argparse.ArgumentParser(
        description='RXNRECer: A deep learning-based tool for predicting enzyme-catalyzed reactions from protein sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic S1 prediction with TSV output
  rxnrecer -i input.fasta -o output.tsv -m s1
  
  # S2 prediction with JSON output
  rxnrecer -i input.fasta -o output.json -m s2 -f json
  
  # S3 prediction with custom batch size
  rxnrecer -i input.fasta -o output.tsv -m s3 -b 50
  
  # Use default output path (temp directory)
  rxnrecer -i input.fasta -m s1
        """
    )
    
    parser.add_argument('-i', '--input_fasta', type=str, default=f'{cfg.DATA_ROOT}sample/sample10.fasta', help='Path to input FASTA file (required)')
    parser.add_argument('-o', '--output_file', type=str, default=f'{cfg.TEMP_DIR}res_sample10.tsv', help='Path to output file (default: temp/res_sample10.tsv)')
    parser.add_argument('-f', '--format', type=str, choices=['tsv', 'json'], default='tsv', help='Output format: tsv or json (default: tsv)')
    parser.add_argument('-m', '--mode', type=str, choices=['s1', 's2', 's3'], default='s1', help='Prediction mode: s1 (basic), s2 (detailed), s3 (LLM reasoning) (default: s1)')
    parser.add_argument('-b', '--batch_size', type=int, default=100, help='Batch size for processing (default: 100)')
    parser.add_argument('-v', '--version', action='version', version='RXNRECer 1.4.0')
    
    # 显示帮助信息
    if len(sys.argv) == 1:
        parser.print_help()
        return
    
    args = parser.parse_args()
    args.output_file = normalize_output_file(args.output_file, args.format)
    
    # 2. 验证输入参数
    if not args.input_fasta:
        print("Error: Input FASTA file is required!")
        print("Use -i or --input_fasta to specify the input file.")
        print("Use -h or --help for more information.")
        return
    
    if not os.path.exists(args.input_fasta):
        print(f"Error: Input file '{args.input_fasta}' does not exist!")
        return
    
    # 3. 准备输出目录
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 4. 显示运行信息
    print("RXNRECer v1.4.0 - Enzyme Reaction Prediction")
    print(f"Input file: {args.input_fasta}")
    print(f"Output file: {args.output_file}")
    print(f"Output format: {args.format}")
    print(f"Prediction mode: {args.mode}")
    print(f"Batch size: {args.batch_size}")
    print("-" * 50)
    
    # 5. 检查缓存
    cache_file_name = None
    if cfg.CACHE_ENABLED:
        cache_file_name = ftool.get_cache_filename(input_file=args.input_fasta, mode=args.mode, output_format=args.format)
        if cache_file_name and ftool.check_cache(cache_file_name):
            res = ftool.load_from_cache(cache_file_name)
            if res is not None:
                save_data(res, args.output_file, args.format)
                print(f"✅ Results loaded from cache and saved to {args.output_file}")
                print(f"⏱️  Total time: {time.perf_counter() - start_time:.2f} seconds")
                return 0
    
    # 6. 验证S3模式配置
    if args.mode == 's3' and (not cfg.LLM_API_KEY or not cfg.LLM_API_URL):
        print("Error: LLM API key and URL are required for S3 mode!")
        return
    
    # 7. 执行预测
    try:
        res = step_by_step_prediction(
            input_data=args.input_fasta, 
            mode=args.mode, 
            batch_size=args.batch_size
        )
        
        # 8. 保存结果
        save_data(res, args.output_file, args.format)
        
        # 9. 保存到缓存
        if cfg.CACHE_ENABLED and cache_file_name:
            ftool.save_to_cache(res, cache_file_name)
        
        # 10. 完成
        elapsed_time = time.perf_counter() - start_time
        print("✅ Prediction completed successfully!")
        print(f"⏱️  Total time: {elapsed_time:.2f} seconds")
        print(f"📁 Results saved to: {args.output_file}")
        
    except Exception as e:
        print(f"❌ Error during prediction: {str(e)}")
        print("Please check your input parameters and try again.")
        return 1
    
    return 0

if __name__ == '__main__':
    main()
    

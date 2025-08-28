import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/')
from rxnrecer.config import config as cfg
import pandas as pd
import numpy as np
import time
from types import SimpleNamespace
import torch
import argparse
from rxnrecer.lib.ml import mlpredict as predRXN
from rxnrecer.lib.llm import qa as llm_qa
from rxnrecer.utils import file_utils as ftool
from rxnrecer.utils import format_utils
from rxnrecer.utils import bio_utils as butils
from tqdm import tqdm
from collections import defaultdict
from rxnrecer.lib.model import mactive as Mactive


def load_model(model_weight_path= cfg.FILE_MOLEL_PRODUCTION_BEST_MODEL):
    mcfg = SimpleNamespace(
    #模型参数
    batch_size=1,
    esm_out_dim=1280,
    gru_h_dim=512,
    att_dim=32,
    dropout_rate=0.2,
    freeze_esm_layers = 32, #冻结层数,
    output_dimensions=10479,
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu"),
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



#region 1. 加载计算数据
def load_data(input_data):
    """
    加载计算数据，支持两种输入格式：
      1) 传入 FASTA 文件路径 (str) -> 自动解析为 DataFrame
      2) 直接传入 DataFrame -> 内含 'seq' 和 'uniprot_id' 列
    """
    if isinstance(input_data, str):
        # 假设传入的是FASTA文件路径
        print(f'Detected input is a FASTA file: {input_data}')
        input_df = ftool.fasta_to_dataframe(fasta_file=input_data)
    elif isinstance(input_data, pd.DataFrame):
        input_df = input_data.copy().reset_index(drop=True)
    else:
        raise ValueError("input_data should be either a path to FASTA or a pandas DataFrame.")
    
    # 统一保证这些列存在（若你的 fasta_to_dataframe 中列名不同，需做相应处理）
    if 'seq' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'seq'.")
    if 'uniprot_id' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'uniprot_id'.")
    
    return input_df
#endregion 


def format_obj(x, ndigits=6):
    """递归处理单元格内容，保留浮点数到指定小数位"""
    if isinstance(x, (np.floating, float)):
        return round(float(x), ndigits)
    elif isinstance(x, dict):
        return {k: format_obj(v, ndigits) for k, v in x.items()}
    elif isinstance(x, (list, tuple)):
        return [format_obj(v, ndigits) for v in x]
    else:
        return x


#region 2. 保存计算结果
def save_data(resdf, output_file, output_format='tsv'):
    """
    保存结果到文件 (TSV/JSON)
    resdf: DataFrame 格式的结果
    output_file: 输出文件路径
    output_format: 'tsv' 或 'json'
    """
    resdf = resdf.applymap(lambda x: format_obj(x, 4))
    if output_format == 'tsv':
        resdf.to_csv(output_file, index=False, sep='\t', float_format="%.4f")
    elif output_format == 'json':
        resdf.to_json(output_file, orient='records', indent=4)
    else:
        print(f'Error: Invalid output format {output_format}. Skip saving.')
#endregion 

def step_by_step_prediction(input_data,
                               output_file=None,
                               output_format='tsv',
                               mode='s1',
                               batch_size=100):
    
    if mode == 's3':
        if not cfg.LLM_API_KEY or not cfg.LLM_API_URL:
            print("Error: LLM API key and URL are required for S3 mode!")
            return
    
    input_df = load_data(input_data)
    print(f'Step 1: Preparing input data, loading {len(input_df)} proteins')

    # Step 2: Load predictive model
    print(f'Step 2: Loading predictive model')
    model, mcfg = load_model(model_weight_path=cfg.FILE_MOLEL_PRODUCTION_BEST_MODEL) 
    
        
    res = []
    
    print(f'Step 3: Running prediction on {len(input_df)} proteins')
    # Total number of batches needed (rounding up for incomplete batches)
    batch_num = (len(input_df) + batch_size - 1) // batch_size
    
    for i in tqdm(range(batch_num)):
        # Calculate batch start and end indices
        start = i * batch_size
        end = min((i + 1) * batch_size, len(input_df))  # Ensure the last batch is handled correctly
        
        # Slice the input_data
        batch_data = input_df.iloc[start:end]
        
        # If the batch is non-empty, run prediction
        if not batch_data.empty:
            res_batch = single_batch_run_prediction(batch_data, model=model, mcfg=mcfg, mode=mode)
            res = res + [res_batch]  # Use extend for better performance
    fres = pd.concat(res, axis=0, ignore_index=True)
    
    if output_file is not None:
        print(f'Step 4: Saving results to {output_file} (format={output_format})')
        save_data(fres, output_file, output_format)
        
    fres = fres.drop_duplicates(subset=['input_id'], keep='first')
    # fres = input_df.rename(columns={'uniprot_id':'input_id'}).merge(fres, on='input_id', how='left')
    return fres
    


#region RXNRECer 预测API
def single_batch_run_prediction(input_df, model, mcfg,  mode='s1'):
    """
    对输入的 DataFrame 进行 RXNRECer 预测，并返回预测结果 DataFrame。
    可选参数:
      - input_df: 输入的 DataFrame，包含 'seq' 和 'uniprot_id' 列
      - model: 预训练好的模型
      - mcfg: 模型参数配置
      - featureBank: 用于计算相似度的特征库
      - getEquation: 是否查询并附加方程信息 (需要 cfg.FILE_DS_RHEA_REACTIONS)
      - Ensemble: 是否使用集成学习方法
    """
    input_df = input_df.reset_index(drop=True)
    # RXNRECer-S1
    print('Running RXNRECer-S1 ...')
    res, res_prob = Mactive.predict_sequences(
        model=model,
        sequences=input_df.seq,
        model_weight_path=mcfg.model_weight_path,
        dict_path=mcfg.dict_path,
        batch_size=2,
        device=mcfg.device
    )
    # 整合预测结果
    # 注意：如果想保留原来所有列，可以在 input_df.copy() 的基础上再拼结果
    res_df_s1 = input_df[['uniprot_id']].reset_index(drop=True).copy()
    res_df_s1['RXNRECer'] = res
    res_df_s1['RXNRECer_with_prob'] = res_prob
    res_df_s1 = res_df_s1.rename(columns={'uniprot_id': 'input_id'})
    
    if mode == 's1':
        #获取反应详情
        rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
        res_df_s1= butils.get_rxn_details_batch(df_rxns=res_df_s1, rxn_bank=rxn_bank, rxn_id_column='RXNRECer')
        res_df_s1['rxn_details'] = res_df_s1.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, RXN_details=x.RXN_details, mode='s2'), axis=1).tolist()
        res_df_s1 = res_df_s1[['input_id','RXNRECer', 'RXNRECer_with_prob', 'rxn_details']]
        return res_df_s1


    # 集成学习
    if mode == 's2':
        print('Running RXNRECer-S2 ...')
        res_df_s2 = get_ensemble(input_df=input_df[['uniprot_id', 'seq']], rxnrecer_df=res_df_s1)#集成学习
        #获取反应详情
        rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)
        res_df_s2= butils.get_rxn_details_batch(df_rxns=res_df_s2, rxn_bank=rxn_bank, rxn_id_column='RXNRECer')
        #格式化输出
        res_df_s2['rxn_details'] = res_df_s2.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, RXN_details=x.RXN_details, mode='s2'), axis=1).tolist()
        res_df_s2 = res_df_s2[['input_id','RXNRECer', 'RXNRECer_with_prob', 'rxn_details']]
        return res_df_s2
    
    if mode == 's3':
        print('Running RXNRECer-S3 ...')
        res_df_s2 = get_ensemble(input_df=input_df[['uniprot_id', 'seq']], rxnrecer_df=res_df_s1).rename(columns={'RXNRECer': 'RXNRECer-S2'})
        s3_input_df = res_df_s2.merge(res_df_s1[['input_id', 'RXNRECer']].rename(columns={'RXNRECer': 'RXNRECer-S1'}), on='input_id', how='left'
                            ).merge(input_df.rename(columns={'uniprot_id': 'input_id'}), on='input_id', how='left')
        res_df_s3 = llm_qa.batch_chat(res_rxnrecer=s3_input_df, api_key=cfg.LLM_API_KEY, api_url=cfg.LLM_API_URL)
        res_df_s3['rxn_details'] = res_df_s3.apply(lambda x: format_utils.format_rxn_output(RXNRECer_with_prob=x.RXNRECer_with_prob, 
                                                                                            RXNRECER_S3=x['RXNRECER-S3'], 
                                                                                            RXN_details=x.RXN_details, mode='s3'), axis=1)
        res_df_s3 = res_df_s3[['input_id','RXNRECer-S2', 'RXNRECer_with_prob', 'rxn_details']].rename(columns={'RXNRECer-S2': 'RXNRECer'})
        
        return res_df_s3
        
#endregion

#region 集成时处理酶和非酶混合情况
def res_refinement(rxn_prob):
    """
    精炼酶预测反应结果，规则：
    1. 只有一个结果（无论是否为非酶），直接保留
    2. 如果 '-' 是最大概率项，则认为是非酶，保留 '-'
    3. 如果 '-' 是最小概率项或居中（非最大非最小），则删除 '-'
    """
    if len(rxn_prob) == 1:
        return cfg.SPLITER.join(rxn_prob.keys()), rxn_prob

    sorted_items = sorted(rxn_prob.items(), key=lambda x: x[1], reverse=True)

    # 如果 '-' 是最高概率项 → 非酶，保留
    if sorted_items[0][0] == '-':
        return '-', {'-': rxn_prob['-']}

    # 如果 '-' 存在且不是最高项 → 无论居中或最小，统一删除
    if '-' in rxn_prob:
        rxn_prob.pop('-')

    return cfg.SPLITER.join(rxn_prob.keys()), rxn_prob
#endregion


def get_ensemble(input_df, rxnrecer_df):
    """
    融合多个方法对蛋白进行反应预测，生成集成预测结果，并处理酶/非酶混合情况。

    输入参数：
    - input_df: 待预测的蛋白信息（DataFrame），需包含 'uniprot_id' 列
    - rxnrecer_df: RXNRECer 预测结果（DataFrame），需包含 'input_id' 和 'RXNRECer_with_prob' 列

    返回：
    - res_df: 仅包含 input_id、融合预测字符串（RXNRECer）和概率字典（RXNRECer_with_prob）的 DataFrame
    """

    # 调用各个预测方法（MSA、CatFam、ECRECer、T5、RXNRECer）
    res_msa = predRXN.getmsa(df_test=input_df, k=1)
    res_catfam = predRXN.getcatfam(df_test=input_df)
    res_ecrecer = predRXN.getecrecer(df_test=input_df)
    res_t5 = predRXN.getT5(df_test=input_df, topk=1)
    res_rxnrecer = rxnrecer_df.copy().rename(columns={'input_id': 'uniprot_id'})

    # 融合不同模型的结果到同一 DataFrame 中（按 uniprot_id 左连接）
    baggingdf = res_rxnrecer.merge(
                    res_msa[['uniprot_id', 'rxn_msa']], on='uniprot_id', how='left'
                ).merge(
                    res_catfam[['uniprot_id', 'rxn_catfam']], on='uniprot_id', how='left'
                ).merge(
                    res_ecrecer[['uniprot_id', 'rxn_ecrecer']], on='uniprot_id', how='left'
                ).merge(
                    res_t5, on='uniprot_id', how='left'
                )

    # 标准化缺失预测：将 NO-PREDICTION、None 和 NaN 统一处理为 '-'
    baggingdf.replace(['NO-PREDICTION', 'None'], '-', inplace=True)
    baggingdf.fillna('-', inplace=True)

    # 调用集成方法进行融合，返回每个蛋白对应的反应预测字典（含概率）
    baggingdf['ensemble'] = baggingdf.apply(lambda x: integrateEnsemble(
        esm=x.esm,
        t5=x.t5,
        rxn_recer=x.RXNRECer_with_prob,
        rxn_msa=x.rxn_msa,
        rxn_catfam=x.rxn_catfam,
        rxn_ecrecer=x.rxn_ecrecer
    ), axis=1)

    # 提取融合结果为两列：
    # - RXNRECer：拼接 reaction ids（字符串）
    # - RXNRECer_with_prob：原始融合字典
    baggingdf['RXNRECer_with_prob'] = baggingdf['ensemble']
    baggingdf['RXNRECer'] = baggingdf['ensemble'].apply(lambda x: cfg.SPLITER.join(x.keys()))

    # 恢复 input_id 命名以便与外部保持一致
    baggingdf.rename(columns={'uniprot_id': 'input_id'}, inplace=True)

    # 仅保留需要输出的字段
    res_df = baggingdf[['input_id', 'RXNRECer', 'RXNRECer_with_prob']].copy()

    #对酶/非酶混合情况进行清洗（去除无效或冲突的 '-' 标签）
    res_df[['RXNRECer', 'RXNRECer_with_prob']] = res_df.apply(
        lambda row: pd.Series(res_refinement(dict(row['RXNRECer_with_prob']))),
        axis=1
    )


    return res_df
    

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
    start_time = time.perf_counter()
    
    # 使用argparse模块解析输入和输出文件路径
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
    
    parser.add_argument('-i', '--input_fasta', 
                       type=str, 
                       help='Path to input FASTA file (required)')
    parser.add_argument('-o', '--output_file', 
                       type=str, 
                       default=f'{cfg.TEMP_DIR}res_sample10.tsv', 
                       help='Path to output file (default: temp/res_sample10.tsv)')
    parser.add_argument('-f', '--format', 
                       type=str, 
                       choices=['tsv', 'json'], 
                       default='tsv', 
                       help='Output format: tsv or json (default: tsv)')
    parser.add_argument('-m', '--mode', 
                       type=str, 
                       choices=['s1', 's2', 's3'], 
                       default='s1', 
                       help='Prediction mode: s1 (basic), s2 (detailed), s3 (LLM reasoning) (default: s1)')
    parser.add_argument('-b', '--batch_size', 
                       type=int, 
                       default=100, 
                       help='Batch size for processing (default: 100)')
    parser.add_argument('-v', '--version', 
                       action='version', 
                       version='RXNRECer 1.0.0')
    
    # 如果没有提供任何参数，显示帮助
    if len(sys.argv) == 1:
        parser.print_help()
        return
    
    args = parser.parse_args()
    
    # 检查必需参数
    if not args.input_fasta:
        print("Error: Input FASTA file is required!")
        print("Use -i or --input_fasta to specify the input file.")
        print("Use -h or --help for more information.")
        return
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input_fasta):
        print(f"Error: Input file '{args.input_fasta}' does not exist!")
        return
    
    # 创建输出目录（如果不存在）
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print(f"RXNRECer v1.0.0 - Enzyme Reaction Prediction")
    print(f"Input file: {args.input_fasta}")
    print(f"Output file: {args.output_file}")
    print(f"Output format: {args.format}")
    print(f"Prediction mode: {args.mode}")
    print(f"Batch size: {args.batch_size}")
    print("-" * 50)
    
    
    if args.mode == 's3':
        if not cfg.LLM_API_KEY or not cfg.LLM_API_URL:
            print("Error: LLM API key and URL are required for S3 mode!")
            return
    
    try:
        res = step_by_step_prediction(
            input_data=args.input_fasta, 
            output_file=args.output_file, 
            output_format=args.format,  
            mode=args.mode, 
            batch_size=args.batch_size
        )
        
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        
        print(f"\n✅ Prediction completed successfully!")
        print(f"⏱️  Total time: {elapsed_time:.2f} seconds")
        print(f"📁 Results saved to: {args.output_file}")
        
    except Exception as e:
        print(f"\n❌ Error during prediction: {str(e)}")
        print("Please check your input parameters and try again.")
        return 1
    
    return 0

if __name__ == '__main__':
    main()
    

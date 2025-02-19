import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')
from config import conf as cfg
import pandas as pd
import numpy as np
from types import SimpleNamespace
from Bio import SeqIO
import torch
import argparse
import hashlib
import modules.simi_caculator as simitool
from modules.predict import predRXN
from tqdm import tqdm
from collections import defaultdict
from methods import Mactive
import json

#region 读取 FASTA 文件并转换为 DataFrame
def fasta_to_dataframe(fasta_file):
    # 使用列表推导式简化数据提取
    data = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    
    # 创建 DataFrame
    df = pd.DataFrame(data, columns=["uniprot_id", "seq"])
    return df
#endregion 


def clean_string(s):
    if isinstance(s, str):
        # 去除多余的转义字符
        return s.replace('\\\\', '\\')
    return s

def clean_data(data):
    if isinstance(data, dict):
        # 递归处理字典中的每一个项
        return {key: clean_data(value) for key, value in data.items()}
    elif isinstance(data, list):
        # 递归处理列表中的每一个项
        return [clean_data(item) for item in data]
    else:
        # 处理单一值
        return clean_string(data)

def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)  # 读取并解析 JSON 文件
        data = clean_data(data)  # 递归清理数据中的转义字符
    return data

def get_rxn_details_from_rxn_json(rxn_ids):
    
    if rxn_ids =='-':
        return '-'
    rxn_ids = rxn_ids.replace(":", "_")  # 去除空格
    rxn_id_array = rxn_ids.split(cfg.SPLITER)
    rxn_list = []  # 用于存储每个DataFrame
    
    for rxn_id in rxn_id_array:
        item = read_json_file(f"{cfg.DIR_PROJECT_ROOT}/{cfg.DIR_RXN_JSON}{rxn_id}.json")
        rxn_list.append(item)
    res = pd.json_normalize(rxn_list)
    return res



def hash_dataframe(df):
    # 将 DataFrame 转换为字符串格式
    df_str = df.to_string(index=False, header=True)
    # 创建哈希对象
    return hashlib.md5(df_str.encode('utf-8')).hexdigest()



def load_model(model_weight_path= cfg.FILE_WEIGHT_PRODUCTION_BEST_MODEL):
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
        print(f'  Detected input is a FASTA file: {input_data}')
        input_df = fasta_to_dataframe(fasta_file=input_data)
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

#region 2. 保存计算结果
def save_data(resdf, output_file, output_format='tsv'):
    """
    保存结果到文件 (TSV/JSON)
    resdf: DataFrame 格式的结果
    output_file: 输出文件路径
    output_format: 'tsv' 或 'json'
    """
    if output_format == 'tsv':
        resdf.to_csv(output_file, index=False, sep='\t')
    elif output_format == 'json':
        # 与原 step_by_step_prediction 同理，可再加 rxn_details 之类的字段
        resdf['rxn_details'] = resdf['RXNRECer'].apply(
            lambda x: get_rxn_details_from_rxn_json(rxn_ids=x) if x != '-' else '-'
        )
        resdf.to_json(output_file, orient='records', indent=4)
    else:
        print(f'Error: Invalid output format {output_format}. Skip saving.')
#endregion 

def step_by_step_prediction(input_data,
                               output_file=None,
                               output_format='tsv',
                               getEquation=False,
                               Ensemble=False,
                               batch_size=100):
    
    print('Step 1: Preparing input data')
    input_df = load_data(input_data)

    # Step 3: Load predictive model
    print('Step 2: Loading predictive model')
    model, mcfg = load_model() 
    
        
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
            res_batch = single_batch_run_prediction(batch_data,
                                                model=model,
                                                mcfg=mcfg,
                                                getEquation=getEquation,
                                                 Ensemble=Ensemble)
            res = res + [res_batch]  # Use extend for better performance
    fres = pd.concat(res, axis=0, ignore_index=True)
    

    
    
    if output_file is not None:
        print(f'Step 5: Saving results to {output_file} (format={output_format})')
        save_data(fres, output_file, output_format)
    return fres
    


#region RXNRECer 预测API
def single_batch_run_prediction(input_df, model, mcfg, getEquation=False, Ensemble=False):
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
    resdf = input_df[['uniprot_id']].reset_index(drop=True).copy()
    resdf['RXNRECer'] = res
    resdf['RXNRECer_with_prob'] = res_prob
    resdf = resdf.rename(columns={'uniprot_id': 'input_id'})
        
        
    # 如果需要从“Rhea reaction 数据库”中查方程式，则执行
    if getEquation:
        print('Fetching reaction equations (RHEA) ...')
        RXNs = pd.read_feather(cfg.FILE_DS_RHEA_REACTIONS)  # 需提前定义 cfg.FILE_DS_RHEA_REACTIONS
        # RXNs 中需至少包含：reaction_id, equation, equation_chebi, ...
        def fetch_equations(rxn_string):
            if rxn_string == '-':
                return '-'
            rxn_list = rxn_string.split(cfg.SPLITER)
            eqs = RXNs[RXNs.reaction_id.isin(rxn_list)].equation.to_list()
            return eqs if eqs else '-'

        def fetch_equations_chebi(rxn_string):
            if rxn_string == '-':
                return '-'
            rxn_list = rxn_string.split(cfg.SPLITER)
            eqs_chebi = RXNs[RXNs.reaction_id.isin(rxn_list)].equation_chebi.to_list()
            return eqs_chebi if eqs_chebi else '-'
        
        resdf['equations'] = resdf['RXNRECer'].apply(fetch_equations)
        resdf['equations_chebi'] = resdf['RXNRECer'].apply(fetch_equations_chebi)
        
        # 集成学习
    if Ensemble:
        resdf = get_ensemble(input_df=input_df[['uniprot_id', 'seq']], rxnrecer_df=resdf)
    return resdf
#endregion




def get_ensemble(input_df, rxnrecer_df):
        res_msa = predRXN.getmsa(df_test=input_df, k=1)
        res_catfam = predRXN.getcatfam(df_test=input_df)
        res_ecrecer = predRXN.getecrecer(df_test=input_df)
        res_rxnrecer = rxnrecer_df.copy().rename(columns={'input_id': 'uniprot_id'})
        res_t5 = predRXN.getT5(df_test=input_df, topk=1)

        baggingdf = res_rxnrecer.merge(
                    res_msa[['uniprot_id', 'rxn_msa']], on='uniprot_id', how='left').merge(
                    res_catfam[['uniprot_id', 'rxn_catfam']], on='uniprot_id', how='left').merge(
                    res_ecrecer[['uniprot_id', 'rxn_ecrecer']], on='uniprot_id', how='left').merge(
                    res_t5, on='uniprot_id', how='left')
        
        baggingdf = baggingdf.replace('NO-PREDICTION', '-')
        baggingdf = baggingdf.replace('None','-')
        baggingdf = baggingdf.fillna('-')
        
        baggingdf['ensemble'] = baggingdf.apply(lambda x: integrateEnsemble(esm=x.esm,
                             t5=x.t5, 
                             rxn_recer=x.RXNRECer_with_prob, 
                             rxn_msa=x.rxn_msa, 
                             rxn_catfam=x.rxn_catfam,
                             rxn_ecrecer=x.rxn_ecrecer), axis=1)
        
       
        baggingdf['RXNRECer']=baggingdf.ensemble.apply(lambda x:(cfg.SPLITER).join(x.keys()))
        baggingdf['RXNRECer_with_prob']=baggingdf.ensemble.copy()
        
        baggingdf = baggingdf.rename(columns={'uniprot_id': 'input_id'})
       
        return baggingdf[['input_id', 'RXNRECer', 'RXNRECer_with_prob']]
    
    

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



        
    

if __name__ == '__main__':
    # start_time = time.perf_counter()
    
    
    # # 使用argparse模块解析输入和输出文件路径
    # parser = argparse.ArgumentParser(description='Run step by step prediction.')
    # parser.add_argument('-i', '--input_fasta', default='data/sample/sample10.fasta', help='Path to input FASTA file')
    # parser.add_argument('-o', '--output_file', default='results/sample/res_sample10.tsv', help='Path to output file')
    # parser.add_argument('-f', '--format', required=False, default='tsv', help='Output format [tsv, json] (default: tsv)')
    # parser.add_argument('-e', '--ensemble', required=False, default=True, help='Whether to use ensemble method (default: False)')
    # parser.add_argument('-g', '--getEquation', required=False, default=False, help='Whether to fetch reaction equations (RHEA) (default: False)')

    # args = parser.parse_args()
    
    # input_fasta = args.input_fasta
    # output_file = args.output_file
    # format = args.format
    # useEsemble = args.ensemble
    # getEquation = args.getEquation

    # res = step_by_step_prediction(input_data=input_fasta, output_file=output_file, output_format=format, getEquation=getEquation, Ensemble=useEsemble,)
    
    # end_time = time.perf_counter()
    # elapsed_time = end_time - start_time
    
    # print(f'Done! \nElapsed time: {end_time - start_time:.2f} seconds')
    
    input_protein_df = pd.read_feather(cfg.FILE_DS_TEST)[['uniprot_id', 'seq']]
    res  = step_by_step_prediction(input_protein_df, cfg.FILE_DS_DICT_RXN2ID, getEquation=False, Ensemble=True, batch_size=1000)
    res.to_pickle(f'{cfg.TEMP_DIR}res_test5years.pkl')
    

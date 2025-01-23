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
from modules.embedding import seqEmbedding as ebdseq
from modules.embedding import t5Embedding as ebdt5
import modules.simi_caculator as simitool
import time
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
def save_data(resdf,output_file, output_format='tsv'):
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
    
    if Ensemble:
        featureBank = pd.read_feather(cfg.FILE_PRODUCTION_FEATURES)
    else:
        featureBank = None
        
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
                                                featureBank=featureBank,
                                                getEquation=getEquation,
                                                Ensemble=Ensemble)
            res = res + [res_batch]  # Use extend for better performance
    fres = pd.concat(res, axis=0, ignore_index=True)
    
    if output_file is not None:
        print(f'Step 5: Saving results to {output_file} (format={output_format})')
        save_data(fres, output_file, output_format)
    return fres
    


#region RXNRECer 预测API
def single_batch_run_prediction(input_df, model, mcfg, featureBank=None, getEquation=False, Ensemble = False):
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
    
    data_hash = hashlib.md5(input_df.to_string().encode('utf-8')).hexdigest()
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer.feather'):
        resdf = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer.feather')
    else:
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
        
        resdf.to_feather(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer.feather')

    # 如果需要从“Rhea reaction 数据库”中查方程式，则执行
    if getEquation:
        print('Step 4.1: Fetching reaction equations (RHEA) ...')
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
        
    if Ensemble:
        path_rxnrecer_ensemble = f'{cfg.TEMP_DIR}{data_hash}_rxnrecer_ensemble.pkl'
        if os.path.exists(path_rxnrecer_ensemble):
            baggingdf = pd.read_pickle(path_rxnrecer_ensemble)
        else:
            res_simi = simi_bagging(input_df, featureBank, topk=1)
            baggingdf = res_simi.merge(resdf, left_on='uniprot_id', right_on='input_id', how='left')
            baggingdf.to_pickle(path_rxnrecer_ensemble)
        
        baggingdf.drop(['input_id'], axis=1, inplace=True)
        baggingdf['ensemble'] = baggingdf.apply(lambda x: integrateEnsemble(x.esm, x.t5, x.RXNRECer_with_prob), axis=1)
        
        return baggingdf

    
    return resdf
#endregion



# Function to calculate similarity between protein features
def get_top_protein_simi(x_feature, y_feature, y_uniprot_id, dict_featureBank,dict_uniprot2rhea, topk):
    future_cosine = simitool.get_cosine_similarity(x_feature, y_feature)        #计算矩阵相似性
    future_cosine = future_cosine.transpose()   
    topk_indices = np.argsort(future_cosine, axis=1)[:, -topk:][:, ::-1]        #查找概率最大的位置
    topk_values = np.round(np.take_along_axis(future_cosine, topk_indices, axis=1), 6)  #找对应位置的概率
    
    #找到对应位置对应的Uniprot_id
    simi_uniprot_id = list(map(dict_featureBank.get, topk_indices.flatten()))
    simi_rhea = list(map(dict_uniprot2rhea.get, simi_uniprot_id))
    result_matrix = np.empty(len(simi_uniprot_id), dtype=object)
    #和概率打包成tuple
    result_matrix[:] = list(zip(simi_rhea, topk_values.flatten()))
    # Reshape to the original shape
    result_matrix = result_matrix.reshape(topk_indices.shape)
    # Create the DataFrame
    res = pd.DataFrame({
        'uniprot_id': y_uniprot_id,
        'simi': list(result_matrix)
    })

    return res


def integrateEnsemble(esm, t5, rxnrecer):
    # Create a dictionary to store all RHEA identifier groups and their probabilities
    merged_dict = defaultdict(list)
    tuples = [esm[0], t5[0]]
    
    # Process each tuple (RHEA IDs and Probabilities)
    for rhea_ids, prob in tuples:
        # Replace None with '-' to ensure they are treated the same
        rhea_ids = rhea_ids if rhea_ids is not None else '-'
        merged_dict[rhea_ids].append(prob)
    
    # Add dictionary data (it already maps RHEA IDs to probabilities)
    for rhea_id, prob in rxnrecer.items():
        # Replace None with '-' for consistency
        rhea_id = rhea_id if rhea_id is not None else '-'
        merged_dict[rhea_id].append(prob)
    
    # For each group of RHEA IDs, we keep the highest probability
    final_result = {}
    for rhea_ids, probs in merged_dict.items():
        final_result[rhea_ids] = max(probs)  # Keep the maximum probability for each group
    
    return final_result


def simi_bagging(input_df, featureBank, topk=3):
    dict_featureBank = pd.Series( featureBank['uniprot_id'],featureBank.index.values).to_dict()
    
    # 从 JSON 文件加载字典数据
    with open(cfg.DICT_UNIPROT_RHEA, "r") as json_file:
        dict_uniprot2rhea = json.load(json_file)
        
    #hash 数据作为文件名    
    data_hash = hashlib.md5(input_df.to_string().encode('utf-8')).hexdigest()
    
    # ESM embedding   
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_esm.feather'):
        embd_esm = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_esm.feather')
    else:
        embd_esm = ebdseq.getEsm(input_df.rename(columns={'uniprot_id':'id'}))
        embd_esm.to_feather(f'{cfg.TEMP_DIR}{data_hash}_esm.feather')
    # ESM similarity
    path_simi_esm = f'{cfg.TEMP_DIR}{data_hash}_simi_esm.pkl'
    if os.path.exists(path_simi_esm):
        esm_cos = pd.read_pickle(path_simi_esm)
    else:
        esm_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.esm), 
                                        y_feature=np.vstack(embd_esm.esm), 
                                        y_uniprot_id=embd_esm.id, 
                                        dict_featureBank=dict_featureBank, 
                                        dict_uniprot2rhea = dict_uniprot2rhea,
                                        topk=topk).rename(columns={'simi':'esm'})
        esm_cos.to_pickle(path_simi_esm)
        
        
    # # Unirep embedding
    # if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather'):
    #     embd_unirep = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather')
    # else:
    #     embd_unirep=ebdseq.getUnirep(input_df.rename(columns={'uniprot_id':'id'}), batch_size=40)
    #     embd_unirep.to_feather(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather')
    # # Unirep similarity
    # path_simi_unirep = f'{cfg.TEMP_DIR}{data_hash}_simi_unirep.pkl'
    # if os.path.exists(path_simi_unirep):
    #     unirep_cos = pd.read_pickle(path_simi_unirep)
    # else:
    #     unirep_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.unirep), 
    #                                     y_feature=np.vstack(embd_unirep.unirep), 
    #                                     y_uniprot_id=embd_unirep.id, 
    #                                     dict_featureBank=dict_featureBank, 
    #                                     dict_uniprot2rhea = dict_uniprot2rhea,
    #                                     topk=topk).rename(columns={'simi':'unirep'})
    #     unirep_cos.to_pickle(path_simi_unirep)
        
        
        
    # T5 Embedding    
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_t5.feather'):
        embd_t5 = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_t5.feather')
    else:
        embd_t5 = ebdt5.get_embd_seq(seqdfwithid=input_df.rename(columns={'uniprot_id':'id'}), batch_szise=20)
        embd_t5.to_feather(f'{cfg.TEMP_DIR}{data_hash}_t5.feather')
    # T5 similarity
    path_simi_t5 = f'{cfg.TEMP_DIR}{data_hash}_simi_t5.pkl'
    if os.path.exists(path_simi_t5):
        t5_cos = pd.read_pickle(path_simi_t5)
    else:
        t5_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.t5), 
                                    y_feature=np.vstack(embd_t5.t5), 
                                    y_uniprot_id=embd_t5.id, 
                                    dict_featureBank=dict_featureBank, 
                                    dict_uniprot2rhea = dict_uniprot2rhea,
                                    topk=topk).rename(columns={'simi':'t5'})
        t5_cos.to_pickle(path_simi_t5)
    
    
    
    res = esm_cos.merge(t5_cos, on='uniprot_id', how='left')
    
    return res
        
    

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
    

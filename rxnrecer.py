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


#region RXNRECer 预测API
def step_by_step_prediction(
    input_data,
    dict_rxn2id,
    output_file=None,
    output_format='tsv',
    getEquation=False,
    Ensemble = False
):
    """
    同时支持:
      1) 传入 FASTA 文件路径 (str) -> 自动解析为 DataFrame
      2) 直接传入 DataFrame -> 内含 'seq' 和 'uniprot_id' 列
    
    可选参数:
      - output_file: 若指定则将结果写入文件 (TSV/JSON)
      - output_format: 'tsv' 或 'json'
      - getEquation: 是否查询并附加方程信息 (需要 cfg.FILE_DS_RHEA_REACTIONS)
    """

    # Step 1: Load reaction encoding dictionary from JSON
    print('Step 1: Load reaction encoding dictionary from JSON file')
    with open(dict_rxn2id, "r") as json_file:
        dict_rxn2id_obj = json.load(json_file)
    print(f'Finished loading rxn2id dictionary. Total {len(dict_rxn2id_obj)} reactions.')

    # Step 2: Load input data -> DataFrame
    print('Step 2: Load input protein data')
    if isinstance(input_data, str):
        # 假设传入的是FASTA文件路径
        print(f'  Detected input is a FASTA file: {input_data}')
        input_df = fasta_to_dataframe(fasta_file=input_data)
    elif isinstance(input_data, pd.DataFrame):
        print('  Detected input is a DataFrame')
        # 为避免 SettingWithCopyWarning，先做一次 copy
        input_df = input_data.copy().reset_index(drop=True)
    else:
        raise ValueError("input_data should be either a path to FASTA or a pandas DataFrame.")

    # 统一保证这些列存在（若你的 fasta_to_dataframe 中列名不同，需做相应处理）
    if 'seq' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'seq'.")
    if 'uniprot_id' not in input_df.columns:
        raise ValueError("Input DataFrame must contain column 'uniprot_id'.")

    
    data_hash = hashlib.md5(input_df.to_string().encode('utf-8')).hexdigest()
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer.feather'):
        resdf = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer.feather')
    else:
        # Step 3: Load predictive model
        print('Step 3: Loading predictive model')
        model, mcfg = load_model() 

        # Step 4: Predict
        print('Step 4: Predicting ...')
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

        
        if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer_ensemble.feather'):
            baggingdf = pd.read_hdf(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer_ensemble.hdf', key='data')
        else:
            res_simi = simi_bagging(input_df , topk=1)
            baggingdf = res_simi.merge(resdf, left_on='uniprot_id', right_on='input_id', how='left')
            baggingdf.to_hdf(f'{cfg.TEMP_DIR}{data_hash}_rxnrecer_ensemble.hdf', key='data', mode='w')
        
        return baggingdf
        baggingdf['ensemble'] = baggingdf.apply(lambda x: integrateEnsemble(x.esm, x.unirep, x.t5, x.RXNRECer_with_prob), axis=1)
        
        return baggingdf

    # 如果需要输出到文件
    if output_file is not None:
        print(f'Step 5: Saving results to {output_file} (format={output_format})')
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


def integrateEnsemble(esm, unirep, t5, rxnrecer):
    # Create a dictionary to store all RHEA identifier groups and their probabilities
    merged_dict = defaultdict(list)
    tuples = [esm[0], unirep[0], t5[0]]
    
    print(esm[0])
    print(unirep[0])
    print(t5[0])
    print(rxnrecer)
    
    
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


def simi_bagging(input_df , topk=3):
    print('  Loading feature Bank')
    featureBank = pd.read_feather(cfg.FILE_PRODUCTION_FEATURES)
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
        print('ESM embedding...')
        embd_esm = ebdseq.getEsm(input_df.rename(columns={'uniprot_id':'id'}))
        embd_esm.to_feather(f'{cfg.TEMP_DIR}{data_hash}_esm.feather')
    # ESM similarity
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_simi_esm.h5'):
        esm_cos = pd.read_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_esm.h5', key='data')
    else:
        print('Caculate ESM simi')
        esm_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.esm), 
                                        y_feature=np.vstack(embd_esm.esm), 
                                        y_uniprot_id=embd_esm.id, 
                                        dict_featureBank=dict_featureBank, 
                                        dict_uniprot2rhea = dict_uniprot2rhea,
                                        topk=topk).rename(columns={'simi':'esm'})
        esm_cos.to_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_esm.h5', key='data', mode='w')
    # Unirep embedding
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather'):
        embd_unirep = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather')
    else:
        print('Unirep embedding')
        embd_unirep=ebdseq.getUnirep(input_df.rename(columns={'uniprot_id':'id'}), batch_size=40)
        embd_unirep.to_feather(f'{cfg.TEMP_DIR}{data_hash}_unirep.feather')
    # Unirep similarity
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_simi_unirep.h5'):
        unirep_cos = pd.read_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_unirep.h5',key='data')
    else:
        print('Caculate Unirep simi')
        unirep_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.unirep), 
                                        y_feature=np.vstack(embd_unirep.unirep), 
                                        y_uniprot_id=embd_unirep.id, 
                                        dict_featureBank=dict_featureBank, 
                                        dict_uniprot2rhea = dict_uniprot2rhea,
                                        topk=topk).rename(columns={'simi':'unirep'})
        unirep_cos.to_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_unirep.h5', key='data', mode='w')
    # T5 Embedding    
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_t5.feather'):
        embd_t5 = pd.read_feather(f'{cfg.TEMP_DIR}{data_hash}_t5.feather')
    else:
        print('T5 embedding...')
        embd_t5 = ebdt5.get_embd_seq(seqdfwithid=input_df.rename(columns={'uniprot_id':'id'}), batch_szise=20)
        embd_t5.to_feather(f'{cfg.TEMP_DIR}{data_hash}_t5.feather')
    # T5 similarity
    if os.path.exists(f'{cfg.TEMP_DIR}{data_hash}_simi_t5.h5'):
        t5_cos = pd.read_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_t5.h5', key='data')
    else:
        print('Caculate T5 simi')
        t5_cos =  get_top_protein_simi(x_feature=np.vstack(featureBank.t5), 
                                    y_feature=np.vstack(embd_t5.t5), 
                                    y_uniprot_id=embd_t5.id, 
                                    dict_featureBank=dict_featureBank, 
                                    dict_uniprot2rhea = dict_uniprot2rhea,
                                    topk=topk).rename(columns={'simi':'t5'})
        t5_cos.to_hdf(f'{cfg.TEMP_DIR}{data_hash}_simi_t5.h5', key='data', mode='w')
    
    
    
    res = esm_cos.merge(unirep_cos, on='uniprot_id', how='left').merge(t5_cos, on='uniprot_id', how='left')
    
    return res
        
    

if __name__ == '__main__':
    start_time = time.perf_counter()
    
    
    # 使用argparse模块解析输入和输出文件路径
    parser = argparse.ArgumentParser(description='Run step by step prediction.')
    parser.add_argument('-i', '--input_fasta', default='data/sample/sample10.fasta', help='Path to input FASTA file')
    parser.add_argument('-o', '--output_file', default='results/sample/res_sample10.tsv', help='Path to output file')
    parser.add_argument('-f', '--format', required=False, default='tsv', help='Output format [tsv, json] (default: tsv)')

    args = parser.parse_args()
    
    input_fasta = args.input_fasta
    output_file = args.output_file
    format = args.format
    
    dict_rxn2id = cfg.FILE_DS_DICT_RXN2ID
   

    res = step_by_step_prediction(input_data=input_fasta, dict_rxn2id = dict_rxn2id, output_file=output_file, output_format=format)
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    
    print(f'Done! \nElapsed time: {end_time - start_time:.2f} seconds')
    

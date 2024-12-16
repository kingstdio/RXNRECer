import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')
from config import conf as cfg
import pandas as pd
from types import SimpleNamespace
from Bio import SeqIO
import torch
import argparse
import hashlib
import time
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
        
        item = read_json_file(f"{cfg.DIR_PROJECT_ROOT}/{cfg.RXN_JSON_DIR}{rxn_id}.json")
        
        rxn_list.append(item)
        
    res = pd.json_normalize(rxn_list)
    
    return res



def hash_dataframe(df):
    # 将 DataFrame 转换为字符串格式
    df_str = df.to_string(index=False, header=True)
    # 创建哈希对象
    return hashlib.md5(df_str.encode('utf-8')).hexdigest()



def load_model():
    
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
    model_weight_path='/hpcfs/fhome/shizhenkun/codebase/preaction/production/files/model_weight/185846best.pth',
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


def step_by_step_prediction(input_fasta, dict_rxn2id,  output_file):
   
    # 从 JSON 文件加载反应编码字典
    print('Step 1: Load reaction encoding dictionary from JSON file')
    
    with open(dict_rxn2id, "r") as json_file:
        dict_rxn2id = json.load(json_file)
        print(f'Finished loading rxn2id dictionary from JSON file. Total {len(dict_rxn2id)} reactions.')  # 打印加载的数据

    #  加载输入反应信息
    print('Step 2: Load input reaction information')
    input_df = fasta_to_dataframe(fasta_file=input_fasta)

    print('Step 3: Loading predictive model')
    model, mcfg = load_model() 

    print('Step 4: Predicting ...')
    res = Mactive.predict_sequences(model=model, 
                            sequences=input_df.seq, 
                            model_weight_path=mcfg.model_weight_path, 
                            dict_path=mcfg.dict_path, 
                            batch_size=2,
                            device=mcfg.device)
    
    input_df['RXNRECer'] = res
    input_df = input_df[['uniprot_id', 'RXNRECer']].rename(columns={'uniprot_id': 'input_id'})
    
    print(f'Step 5: Saving results to {output_file}')
    if format == 'tsv':
        input_df.to_csv(output_file, index=False, sep='\t')
        
    elif format == 'json':
        input_df['rxn_details']=input_df.RXNRECer.apply(lambda x: get_rxn_details_from_rxn_json(rxn_ids=x))
        input_df.to_json(output_file,orient='records', indent=4)
    else:
        print(f'Error: Invalid output format {format}.')
    
    return input_df

#region RXNRECer 预测API
def step_by_step_prediction_with_protein_df(input_protein_df, dict_rxn2id):
   
    # 从 JSON 文件加载反应编码字典
    with open(dict_rxn2id, "r") as json_file:
        dict_rxn2id = json.load(json_file)
        
    model, mcfg = load_model() 
    res = Mactive.predict_sequences(model=model, 
                            sequences=input_protein_df.seq, 
                            model_weight_path=mcfg.model_weight_path, 
                            dict_path=mcfg.dict_path, 
                            batch_size=2,
                            device=mcfg.device)
    
    input_protein_df['RXNRECer'] = res
    input_protein_df = input_protein_df[['uniprot_id', 'RXNRECer']].rename(columns={'uniprot_id': 'input_id'})

    return input_protein_df
#endregion

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
   

    res = step_by_step_prediction(input_fasta=input_fasta, dict_rxn2id = dict_rxn2id, output_file=output_file)
    
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    
    print(f'Done! \nElapsed time: {end_time - start_time:.2f} seconds')
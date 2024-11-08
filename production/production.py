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
from  modules.ml import Mactive
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

def get_rxn_details_from_rxn_json(rxn_ids, rxn_info_base):
    
    rxn_ids = rxn_ids.replace("-;", "").replace(";-", "")  # 去除空格
    if rxn_ids =='-':
        return '-'
    rxn_ids = rxn_ids.replace(":", "_")  # 去除空格
    rxn_id_array = rxn_ids.split(cfg.SPLITER)
    rxn_list = []  # 用于存储每个DataFrame
    
    for rxn_id in rxn_id_array:
        # print(f"{rxn_info_base}{rxn_id}.json")
        item = read_json_file(f"{rxn_info_base}{rxn_id}.json")
        
        rxn_list.append(item)
        
    res = pd.json_normalize(rxn_list)
    
    return res

def extract_equations(cellitem):
    """Extract equations from a cell item"""
    if str(cellitem)=='-':
        return ['-','-']
    
    str_equations = cellitem.reaction_equation.to_list()
    str_equations_chebi = cellitem.reaction_equation_ref_chebi.to_list()
    
    return [str_equations, str_equations_chebi]



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
    model_weight_path= cfg.FILE_MOLEL_PRODUCTION_BEST_MODEL,
    dict_path = cfg.FILE_DS_DICT_ID2RXN
    )


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


def step_by_step_prediction(input_fasta, dict_rxn2id, rxn_info_base, output_file, format='tsv'):
   
    
   
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
    res, res_prob = Mactive.predict_sequences(model=model, 
                            sequences=input_df.seq, 
                            model_weight_path=mcfg.model_weight_path, 
                            dict_path=mcfg.dict_path, 
                            batch_size=2,
                            device=mcfg.device)
    
    # print(res_prob)
    
    input_df['RXNRECer'] = res
    input_df['RXNRECer_with_prob'] = res_prob
    input_df = input_df[['uniprot_id', 'RXNRECer', 'RXNRECer_with_prob']].rename(columns={'uniprot_id': 'input_id'})
    
    print(f'Step 5: Saving results to {output_file}')
    
    # 获取反应详细信息
    input_df['rxn_details']=input_df.RXNRECer.apply(lambda x: get_rxn_details_from_rxn_json(rxn_ids=x, rxn_info_base=rxn_info_base))
    
    if format == 'tsv':
        input_df[['equations', 'equations_chebi']] = input_df.apply(lambda x: extract_equations(x.rxn_details), axis=1, result_type='expand') # 提取反应方程式
        input_df=input_df[['input_id','RXNRECer','RXNRECer_with_prob','equations','equations_chebi']]
        input_df.to_csv(output_file, index=False, sep='\t')
        
    elif format == 'json':
        # input_df['rxn_details']=input_df.RXNRECer.apply(lambda x: get_rxn_details_from_rxn_json(rxn_ids=x))
        input_df.to_json(output_file,orient='records', indent=4)
    else:
        print(f'Error: Invalid output format {format}.')
    
    return input_df


if __name__ == '__main__':
    
    # 使用argparse模块解析输入和输出文件路径
    parser = argparse.ArgumentParser(description='Run step by step prediction.')
    parser.add_argument('-i', '--input_fasta', required=True, help='Path to input FASTA file')
    parser.add_argument('-o', '--output_file', required=True, help='Path to output file')
    parser.add_argument('-f', '--format', required=False, default='tsv', help='Output format [tsv, json] (default: tsv)')

    args = parser.parse_args()
    
    input_fasta = args.input_fasta
    output_file = args.output_file
    format = args.format
    
    dict_rxn2id = cfg.FILE_DS_DICT_RXN2ID
    rxn_info_base = cfg.DIR_RXN_JSON

    res = step_by_step_prediction(input_fasta=input_fasta, 
                                  dict_rxn2id = dict_rxn2id, 
                                  rxn_info_base= rxn_info_base, 
                                  output_file=output_file,
                                  format=format
                                  )
    
    print('Done!')
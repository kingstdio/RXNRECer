'''
Author: Zhenkun Shi
Date: 2023-04-20 06:23:40
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-06-25 16:59:57
FilePath: /RXNRECer/tools/btools.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''

import json
import pandas as pd
import numpy as np
import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')
import requests
from config import conf as cfg
from  tools import evaluation as eva
from tkinter import _flatten



def rxn_eva_metric(eva_df, eva_name, methods, average_type='weighted'):
    print(f'Evaluating: Reaction Predcition Results {eva_name}')
    resl = []
    for m in methods:
        res_item  = eva.caculateMetrix(groundtruth=np.stack(eva_df.lb_rxn_groundtruth), predict=np.stack(eva_df[f'lb_rxn_{m}']), baselineName=m, type='multi', print_flag=False, averege_type=average_type)
        resl.append(res_item)

    resl = pd.DataFrame(resl, columns=['baselineName','mAccuracy','mPrecision','mRecall','mF1'])
    return resl

def rxn_eva_metric_with_colName(eva_df, col_groundtruth, col_pred, eva_name='', average_type='weighted'):
    # print(f'Evaluating: Reaction Predcition Results {eva_name}')
    
    groundtruth = np.stack(eva_df[col_groundtruth])
    pred = np.stack(eva_df[col_pred])
    

    
    res_item  = eva.caculateMetrix(groundtruth=groundtruth, 
                                   predict=pred, 
                                   baselineName=eva_name, 
                                   type='multi', 
                                   print_flag=False, 
                                   averege_type=average_type
                                   )
    # print(res_item)
    # resl = pd.DataFrame(res_item)
    resl = pd.DataFrame([res_item], columns=['mAccuracy','mPrecision','mRecall','mF1', 'avgType'])

    return resl


# 将EC转化为反应
def retrival_reaction_from_ec(ec_pred, ec_reaction_map):
    reaction = '-'
    
    # 非酶无反应
    if ec_pred == '-':
        return '-'
    
    if ec_pred =='NO-PREDICTION':
        return 'NO-PREDICTION'
    
    if ';' not in ec_pred: #对应单个EC
        reaction = ec_reaction_map[ec_reaction_map.ec==ec_pred].reaction_id.to_list() #获取EC对应的反应
        if len(reaction)>0: #如果有对应反应
            reaction = reaction[0].split(cfg.SPLITER) #按分号拆开
    else:   #对应多个EC
        ecs = ec_pred.split(cfg.SPLITER)
        ecs = [item.strip() for item in ecs]

        reaction_temp=[]
        for itemec in ecs:
            reaction_temp  = reaction_temp + ec_reaction_map[ec_reaction_map.ec==itemec].reaction_id.to_list()
        reaction= list(set(_flatten([item.split(cfg.SPLITER) for item in reaction_temp])))
    
    res = (cfg.SPLITER).join(list(set(reaction)))
    
    if res!='':
        return res
    else:
        return 'EC-WITHOUT-REACTION'


#region labelmaker
def make_label(reaction_id, rxn_label_dict):
    """
    将反应id转化成one-hot编码标签，并返回
    输入：
        reaction_id: 反应id，以分号分隔
        rxn_label_dict: 反应id到标签的映射字典
    输出：
        resArray: one-hot编码标签
    """
    # 初始化一个零数组
    resArray = np.zeros(len(rxn_label_dict), dtype=int)
    
    if reaction_id in {'EC-WITHOUT-REACTION', 'NO-PREDICTION'}:  # 如果没有反应或没有预测，随机返回一个标签
        resArray[np.random.randint(len(resArray))] = 1
        return resArray
    
    # 分割反应id并迭代
    for item in reaction_id.split(cfg.SPLITER):
        array_index = rxn_label_dict.get(item)
        if array_index is not None:
            if not isinstance(array_index, int):
                print(f"[WARNING] Invalid index type: {type(array_index)} | value: {array_index} | item: {item}")
            resArray[array_index] = 1

    return resArray

#endregion



#format EC resluts from clean
def load_clean_resluts(res_file):
    data = pd.read_csv(res_file, sep='\t').rename(columns={'Entry':'uniprot_id'})
    def format_ec(eclist):
        res = [item.split('/')[0].replace('EC:','')  for item in eclist.split(';') if item!=None]
        res = cfg.SPLITER.join(res)
        return res
    

    data['ec_clean'] = data.apply(lambda x : format_ec(eclist=x.clean_pred_ec_maxsep), axis=1)
    return data[['uniprot_id', 'ec_clean']]

# region load different method results

def load_deepec_resluts(filepath):
    """load deepec predicted resluts
    Args:
        filepath (string): deepec predicted file
    Returns:
        DataFrame: columns=['id', 'ec_deepec']
    """
    res_deepec = pd.read_csv(f'{filepath}', sep='\t',names=['id', 'ec_number'], header=0 )
    res_deepec.ec_number=res_deepec.apply(lambda x: x['ec_number'].replace('EC:',''), axis=1)
    res_deepec.columns = ['id','ec_deepec']
    res = []
    for index, group in  res_deepec.groupby('id'):
        if len(group)==1:
            res = res + [[group.id.values[0], group.ec_deepec.values[0]]]
        else:
            ecs_str = ';'.join(group.ec_deepec.values)
            res = res +[[group.id.values[0],ecs_str]] 
    res_deepec = pd.DataFrame(res, columns=['id', 'ec_deepec'])
    return res_deepec


def load_praim_res(resfile):
    """[加载PRIAM的预测结果]
    Args:
        resfile ([string]): [结果文件]
    Returns:
        [DataFrame]: [结果]
    """
    f = open(resfile)
    line = f.readline()
    counter =0
    reslist=[]
    lstr =''
    subec=[]
    while line:
        if '>' in line:
            if counter !=0:
                reslist +=[[lstr, ';'.join(subec)]]
                subec=[]
            lstr = line.replace('>', '').replace('\n', '')
        elif line.strip()!='':
            ecarray = line.split('\t')
            subec += [(ecarray[0].replace('#', '').replace('\n', '').replace(' ', '') )]

        line = f.readline()
        counter +=1
    f.close()
    res_priam=pd.DataFrame(reslist, columns=['id', 'ec_priam'])
    return res_priam


def load_catfam_res(resfile):
    res_catfam = pd.read_csv(resfile, sep='\t', names=['id', 'ec_catfam'])
    res_catfam = res_catfam.groupby('id', as_index=False)['ec_catfam'].agg(lambda x: ';'.join(x.dropna())).replace('', '-')
    return res_catfam


def load_ecpred_res(resfile):
    res_ecpred = pd.read_csv(f'{resfile}', sep='\t', header=0)
    res_ecpred = res_ecpred.rename(columns={'Protein ID':'id','EC Number':'ec_ecpred','Confidence Score(max 1.0)':'pident_ecpred'})
    res_ecpred['ec_ecpred']= res_ecpred.ec_ecpred.apply(lambda x : '-' if x=='non Enzyme' else x) 
    return res_ecpred

#endregion




def read_h5_file(file_path):
    with pd.HDFStore(file_path, 'r') as h5:
        data = h5['data']
    return data


def get_simi_Pred(pred_list, uniprot_rxn_dict, topk=3):
    uniprot_id_list = [item[0] for item in pred_list][:topk]
    rxn_ids = [uniprot_rxn_dict.get(uniprot_id) for uniprot_id in uniprot_id_list]
    rxn_res = (cfg.SPLITER).join(set(rxn_ids))
    return rxn_res


def load_dict_rxn2ec():
    
    with open(cfg.DICT_RHEA_EC, 'r') as f:
        dict_rxn2ec = json.load(f)
    return dict_rxn2ec

def load_dict_ec2rxn():
    with open(cfg.DICT_EC_RHEA, 'r') as f:
        dict_ec2rxn = json.load(f)
    return dict_ec2rxn




def transRXN2EC(rxns, dict_rxn2ec):
    if rxns=='-':
        return '-'
    rxn_list = rxns.split(';')
    ec_list = []
    for rxn in rxn_list:
        if rxn in dict_rxn2ec:
            ec_list.append(dict_rxn2ec.get(rxn))
        else:
            ec_list.append('-')
            print(f'{rxn} not in dict_rxn2ec')
    res = ')('.join(ec_list)
    res =f'({res})'
    return res


# 转换为DataFrame
def json_to_dataframe(json_data):
    # 提取UniProt ID
    uniprot_id = list(json_data.keys())[0]
    
    # 为每个条目添加UniProt ID
    for entry in json_data[uniprot_id]:
        entry['uniprot_id'] = uniprot_id
    
    # 创建DataFrame
    df = pd.DataFrame(json_data[uniprot_id])
    
    # 重新排列列顺序，将uniprot_id放在前面
    cols = ['uniprot_id'] + [col for col in df.columns if col != 'uniprot_id']
    df = df[cols]
    
    return df


# 假设df是您提供的DataFrame
def select_best_pdb(df):
    # 规则1+2：先按分辨率升序，再按实验方法排序（X-ray优先）
    df_sorted = df.sort_values(
        by=['resolution', 'experimental_method'],
        ascending=[True, False]  # resolution越小越好，method按字母倒序X-ray优先
    )
    
    # 规则3：如果分辨率相同，选择链ID字母序靠前的（A链优先）
    df_sorted = df_sorted.sort_values('chain_id', ascending=True)
    
    # 规则4（可选）：如果需要特定物种，可以添加筛选
    # df_sorted = df_sorted[df_sorted['tax_id'] == 特定物种ID]
    
    # 返回第一个（最优）条目
    best_row = df_sorted.iloc[0]
    return best_row['pdb_id'], best_row['chain_id']


def download_pdb(pdb_id, save_path=None):
    """下载PDB文件"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        if save_path:
            with open(save_path, 'w') as f:
                f.write(response.text)
            # print(f"PDB文件已保存到: {save_path}")
            return save_path
        return response.text
    except Exception as e:
        print(f"下载失败: {str(e)}")
        return None
    

def parse_pdb_info(pdb_json):
    """解析PDB JSON数据并提取关键信息"""
    # 主要信息提取
    main_info = {
        "pdb_id": pdb_json.get("entry", {}).get("id", ""),
        "resolution": pdb_json.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
        "method": pdb_json.get("exptl", [{}])[0].get("method", "").title(),
        "journal": f"{pdb_json.get('citation', [{}])[0].get('journal_abbrev', '')} "
                  f"({pdb_json.get('citation', [{}])[0].get('year', '')})",
        "ref_doi": pdb_json.get("citation", [{}])[0].get("pdbx_database_id_doi", "")
    }

    # 配体信息提取
    ligands = []
    if "rcsb_binding_affinity" in pdb_json:
        ligands = [{
            "chemical_id": lig["chemical_id"],
            "chemical_name": lig.get("chemical_name", ""),
            "formula": lig.get("formula", "")
        } for lig in pdb_json["rcsb_binding_affinity"]]

    return main_info, ligands

def get_best_pdb(uniprot_id: str):
    """
    修正版：通过 UniProt ID 获取最优 PDB 结构
    返回: {
        "pdb_id": str,
        "resolution": float,
        "method": str,
        "ligands": list
    }
    """
    # 步骤1：通过 PDBe API 获取映射数据
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    response = requests.get(url).json()
    if not response:
        # print (f"No PDB found Use Alphfold2 predicted PDB instead for UniProt ID: {uniprot_id}")
        pdb_path = f'/hpcfs/fpublic/database/alphafold/predicted_pdbs/AF-{uniprot_id}-F1-model_v4.pdb.gz'
        return pdb_path
    
    response = json_to_dataframe(response)
    
    # 执行筛选
    best_pdb, best_chain = select_best_pdb(response)
    # print(f"Best PDB ID: {best_pdb}")
    

    # 步骤4：获取每个 PDB 的元数据

    try:
        pdb_api = f"https://data.rcsb.org/rest/v1/core/entry/{best_pdb}"
        pdb_info = requests.get(pdb_api).json()
        main_info, ligands = parse_pdb_info(pdb_info)
        print(main_info, ligands)
        
    except Exception as e:
        print(f"Error processing PDB RCSB_{uniprot_id}_{best_pdb}: {str(e)}")
    
    pdb_path = f'{cfg.DIR_PDB_BEST}RCSB_UniProt_{uniprot_id}_{best_pdb}.pdb'
    download_pdb(best_pdb, pdb_path)
    return pdb_path
    
    
def get_rxn_detail_by_ids(rxn_id):
    if rxn_id == '-':
        return {
            'reaction id': '-',
            'reaction equation': '-'
        }
    
    rxn_bank =  pd.read_feather(cfg.FILE_RHEA_REACTION)
    rxn_record = rxn_bank[rxn_bank.reaction_id == rxn_id]
    
    if rxn_record.empty:
        return {}  # Avoid errors if the record is missing

    rxn_record = rxn_record.iloc[0]  # Take the first match

    return {
        'reaction id': rxn_record.reaction_id,
        'reaction equation': rxn_record.equation,
        'reaction equation in ChEBI format': rxn_record.equation_chebi,
        'reaction equation in SMILES format': rxn_record.equation_smiles,
        'reaction associated Enzyme Commission Number': rxn_record.ec_number
    }
    
    
def get_rxn_equ_by_ids(rxn_ids, rxn_bank=None):
    if not isinstance(rxn_ids, str) or rxn_ids.strip() == '-':
        return '-'

    rxn_ids = rxn_ids.strip()
    
    if rxn_bank is None:
        rxn_bank = pd.read_feather(cfg.FILE_RHEA_REACTION)

    rxn_list = [rxn_ids] if cfg.SPLITER not in rxn_ids else rxn_ids.split(cfg.SPLITER)
    rxn_record = rxn_bank[rxn_bank.reaction_id.isin(rxn_list)]

    if rxn_record.empty:
        return '-'

    return cfg.SPLITER.join(rxn_record.equation.tolist())



if __name__ =='__main__':
    print('success')
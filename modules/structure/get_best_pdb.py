
import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../../../')
from config import conf as cfg
import pandas as pd
import subprocess
import shutil
from tools import filetool 
from itertools import combinations
from tqdm import tqdm
from Bio.PDB import PDBList
import requests



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
            print(f"PDB文件已保存到: {save_path}")
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

def get_best_pdb(uniprot_id, target_dir=''):
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
        raise ValueError(f"No PDB found for UniProt ID: {uniprot_id}")
    
    response = json_to_dataframe(response)
    
    # 执行筛选
    best_pdb, best_chain = select_best_pdb(response)
    print(f"Best PDB ID: {uniprot_id}->{best_pdb}")
    

    # 步骤4：获取每个 PDB 的元数据

    try:
        pdb_api = f"https://data.rcsb.org/rest/v1/core/entry/{best_pdb}"
        pdb_info = requests.get(pdb_api).json()
        main_info, ligands = parse_pdb_info(pdb_info)
        print(main_info, ligands)
        
    except Exception as e:
        print(f"Error processing PDB RCSB_{uniprot_id}_{best_pdb}: {str(e)}")
       
    download_pdb(best_pdb, f"{target_dir}/RCSB_{uniprot_id}_{best_pdb}.pdb")


    
    # return best_pdb
    

if __name__ == '__main__':
    uniprot_id = "X5F427"
    get_best_pdb(uniprot_id, target_dir='./')
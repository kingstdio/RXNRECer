"""
Biological utility functions for RXNRECer
"""

import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/../')
from rxnrecer.config import config as cfg
from rxnrecer.utils import file_utils as ftool
import re
import subprocess
import tempfile
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML
from pandarallel import pandarallel
# from .file_utils import dataframe_to_fasta  # centralize I/O helpers
import requests
from requests.exceptions import RequestException


def split_compound_str(compounds_str: str, remove_reaction_coefficient: bool = False) -> List[str]:
    """
    Split compound string into list of compounds.
    
    Args:
        compounds_str: String containing compounds
        remove_reaction_coefficient: Whether to remove reaction coefficients
        
    Returns:
        List of compound strings
    """
    if not compounds_str or compounds_str == '-':
        return []
    
    # Split by common separators
    compounds = re.split(r'[+\s]+', compounds_str.strip())
    
    if remove_reaction_coefficient:
        # Remove numeric coefficients
        compounds = [re.sub(r'^\d+\.?\d*\s*', '', comp) for comp in compounds]
    
    # Filter out empty strings
    compounds = [comp.strip() for comp in compounds if comp.strip()]
    
    return compounds


def get_blast_results(train_df: pd.DataFrame, test_df: pd.DataFrame, 
                     k: int = 1, evalue: float = 1e-5) -> pd.DataFrame:
    """
    Perform BLAST search between train and test datasets.
    
    Args:
        train_df: Training dataset DataFrame
        test_df: Test dataset DataFrame
        k: Number of top hits to return
        evalue: E-value threshold
        
    Returns:
        DataFrame with BLAST results
    """
    with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_train, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.tsv') as res_blast, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.dmnd') as db_dmnd:
        
        # Convert DataFrames to FASTA
        dataframe_to_fasta(train_df, fasta_train.name)
        dataframe_to_fasta(test_df, fasta_test.name)
        
        # Build DIAMOND database
        cmd1 = ["diamond", "makedb", "--in", fasta_train.name, "-d", db_dmnd.name, "--quiet"]
        subprocess.run(cmd1, check=True)
        
        # Run BLAST search
        cmd2 = ["diamond", "blastp", "-d", db_dmnd.name, "-q", fasta_test.name, 
                "-o", res_blast.name, "-b5", "-c1", "-k", str(k), 
                "-e", str(evalue), "--quiet"]
        subprocess.run(cmd2, check=True)
        
        # Read results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        results = pd.read_csv(res_blast.name, sep='\t', names=columns)
        
    return results


# Removed duplicated dataframe_to_fasta to avoid divergence; use file_utils.dataframe_to_fasta


def calculate_sequence_similarity(seq1: str, seq2: str, method: str = 'identity') -> float:
    """
    Calculate sequence similarity between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        method: Similarity method ('identity', 'blosum62', 'pam250')
        
    Returns:
        Similarity score
    """
    if method == 'identity':
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)
    else:
        # TODO: Implement other similarity methods
        raise NotImplementedError(f"Method {method} not implemented")


def get_protein_properties(sequence: str) -> Dict[str, float]:
    """
    Calculate basic protein properties.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Dictionary with protein properties
    """
    # Amino acid frequencies
    aa_freq = {}
    for aa in sequence:
        aa_freq[aa] = aa_freq.get(aa, 0) + 1
    
    # Normalize frequencies
    seq_len = len(sequence)
    aa_freq = {aa: count/seq_len for aa, count in aa_freq.items()}
    
    # Calculate properties
    properties = {
        'length': seq_len,
        'molecular_weight': calculate_molecular_weight(sequence),
        'isoelectric_point': calculate_isoelectric_point(sequence),
        'hydrophobicity': calculate_hydrophobicity(sequence),
        'charge': calculate_charge(sequence)
    }
    
    return properties


def calculate_molecular_weight(sequence: str) -> float:
    """Calculate molecular weight of protein sequence."""
    # Amino acid molecular weights (Da)
    aa_weights = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    
    weight = sum(aa_weights.get(aa, 0) for aa in sequence.upper())
    return weight - 18.0 * (len(sequence) - 1)  # Subtract water molecules


def calculate_isoelectric_point(sequence: str) -> float:
    """Calculate isoelectric point of protein sequence."""
    # pKa values for amino acids
    pka_values = {
        'D': 3.65, 'E': 4.25, 'H': 6.0, 'K': 10.53, 'R': 12.48,
        'Y': 10.07, 'C': 8.18
    }
    
    # Count charged residues
    charged_residues = {aa: sequence.upper().count(aa) for aa in pka_values}
    
    # Simple calculation (can be improved with more sophisticated methods)
    net_charge = (charged_residues['D'] + charged_residues['E'] - 
                  charged_residues['K'] - charged_residues['R'] - 
                  charged_residues['H'])
    
    # Approximate pI calculation
    if net_charge > 0:
        return 10.0
    elif net_charge < 0:
        return 4.0
    else:
        return 7.0


def calculate_hydrophobicity(sequence: str) -> float:
    """Calculate hydrophobicity of protein sequence."""
    # Kyte-Doolittle hydrophobicity scale
    hydrophobicity_scores = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    scores = [hydrophobicity_scores.get(aa, 0) for aa in sequence.upper()]
    return sum(scores) / len(scores)


def calculate_charge(sequence: str) -> float:
    """Calculate net charge of protein sequence at pH 7."""
    # Charge at pH 7
    charges = {
        'D': -1, 'E': -1, 'H': 0.1, 'K': 1, 'R': 1, 'Y': 0
    }
    
    net_charge = sum(charges.get(aa, 0) for aa in sequence.upper())
    return net_charge


def download_uniprot_entry(uniprot_id: str, save_path: Optional[str] = None) -> str:
    """
    Download UniProt entry.
    
    Args:
        uniprot_id: UniProt ID
        save_path: Path to save file (optional)
        
    Returns:
        Path to downloaded file
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    
    if save_path is None:
        save_path = f"{uniprot_id}.fasta"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        with open(save_path, 'w') as f:
            f.write(response.text)
        
        return save_path
    except RequestException as e:
        raise Exception(f"Failed to download UniProt entry {uniprot_id}: {e}")


def parse_fasta_header(header: str) -> Dict[str, str]:
    """
    Parse FASTA header to extract information.
    
    Args:
        header: FASTA header line
        
    Returns:
        Dictionary with parsed information
    """
    info = {}
    
    # Extract ID
    if '|' in header:
        parts = header.split('|')
        if len(parts) >= 2:
            info['database'] = parts[0].lstrip('>')
            info['id'] = parts[1]
            if len(parts) >= 3:
                info['entry_name'] = parts[2]
    
    # Extract description
    if ' ' in header:
        description = header.split(' ', 1)[1]
        info['description'] = description
    
    return info


def validate_sequence(sequence: str) -> bool:
    """
    Validate protein sequence.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        True if valid, False otherwise
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in sequence.upper())


def clean_sequence(sequence: str) -> str:
    """
    Clean protein sequence by removing invalid characters.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Cleaned sequence
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return ''.join(aa for aa in sequence.upper() if aa in valid_aa)

#region 根据反应ID获取反应详细信息
def get_rxn_detail(rxn_id, rxn_bank):
    """
    根据反应ID获取单个反应的详细信息，返回Reaction对象
    
    Args:
        rxn_id (str): 反应ID字符串，如 'RHEA:12345'
        rxn_bank (pd.DataFrame): 反应数据库DataFrame，需包含列：
            - reaction_id: 反应ID
            - equation: 反应方程式
            - equation_chebi: ChEBI格式的反应方程式
            - equation_smiles: SMILES格式的反应方程式
            - ec_number: 酶委员会编号
        
    Returns:
        Reaction: Reaction对象，包含完整的反应信息
            如果反应ID无效或不存在，返回None
            
    Note:
        使用Reaction对象可以统一接口，方便后续处理
    """
    from rxnrecer.lib.rxn.Reaction import Reaction
    
    if not rxn_id or rxn_id == '-':
        return None

    match = rxn_bank[rxn_bank.reaction_id == rxn_id.strip()]
    if match.empty:
        return None

    row = match.iloc[0]
    
    try:
        # 创建Reaction对象
        reaction = Reaction(
            rxn_smiles=row.equation_smiles,
            rxn_equation=row.equation,
            rxn_equation_ref_chebi=row.equation_chebi,
            rxn_id=row.reaction_id,
            rxn_ec=row.ec_number
        )
        return reaction
    except Exception as e:
        print(f"警告: 无法创建反应对象 {rxn_id}: {e}")
        return None
#endregion

def get_rxn_details_from_rxn_json(rxn_ids):
    """
    从JSON文件获取反应详情，支持多个分隔符和异常处理
    
    Args:
        rxn_ids (str): 反应ID字符串，支持多种分隔符（; | ,）
            例如: 'RHEA:14709;RHEA:24076;RHEA:32187' 或 'RHEA:14709|RHEA:24076'
        
    Returns:
        pd.DataFrame: 包含反应详情的DataFrame，如果出错返回空DataFrame
        
    Note:
        自动处理多种分隔符，支持 ; | , 等分隔符
        对找不到的JSON文件进行异常处理，跳过无效文件
    """
    import pandas as pd
    
    if not rxn_ids or rxn_ids == '-':
        return pd.DataFrame()
    
    # 处理多种分隔符，支持 ; | , 等
    separators = [';', '|', ',', cfg.SPLITER]
    rxn_id_array = []
    
    # 尝试不同的分隔符分割
    for sep in separators:
        if sep in rxn_ids:
            rxn_id_array = [rxn_id.strip() for rxn_id in rxn_ids.split(sep) if rxn_id.strip()]
            break
    
    # 如果没有找到分隔符，当作单个ID处理
    if not rxn_id_array:
        rxn_id_array = [rxn_ids.strip()]
    
    rxn_list = []  # 用于存储每个JSON数据
    
    for rxn_id in rxn_id_array:
        try:
            # 处理RHEA:格式的ID，转换为文件名格式
            if ':' in rxn_id:
                file_id = rxn_id.replace(":", "_")
            else:
                file_id = rxn_id
            
            # 构建文件路径
            file_path = f"{cfg.DIR_RXN_JSON}{file_id}.json"
            
            # 检查文件是否存在
            if not os.path.exists(file_path):
                print(f"警告: 找不到反应文件 {file_path}")
                continue
            
            # 读取JSON文件
            item = ftool.read_json_file(file_path)
            if item:  # 确保读取成功
                rxn_list.append(item)
            else:
                print(f"警告: 无法读取反应文件 {file_path}")
                
        except Exception as e:
            print(f"警告: 处理反应ID {rxn_id} 时出错: {e}")
            continue
    
    # 如果没有成功读取任何文件，返回空DataFrame
    if not rxn_list:
        print(f"警告: 没有成功读取任何反应文件")
        return pd.DataFrame()
    
    try:
        # 使用pandas json_normalize处理数据
        res = pd.json_normalize(rxn_list)
        return res
    except Exception as e:
        print(f"警告: 数据标准化失败: {e}")
        return pd.DataFrame()


def get_rxn_details_list(rxn_string, rxn_bank, spliter=cfg.SPLITER):
    """
    解析包含多个反应ID的字符串，返回每个反应的详细信息列表
    
    Args:
        rxn_string (str): 反应字符串，可能包含多个反应ID，用分隔符分隔
            例如: 'RHEA:12345|RHEA:67890' 或 'RHEA:12345'
        rxn_bank (pd.DataFrame): 反应数据库DataFrame
        spliter (str, optional): 分隔符，默认为配置文件中的SPLITER
        
    Returns:
        list: Reaction对象的列表，每个对象包含完整的反应信息
        
    Examples:
        >>> get_rxn_details_list('RHEA:12345|RHEA:67890', rxn_bank, '|')
        [<Reaction object>, <Reaction object>]
        
        >>> get_rxn_details_list('-', rxn_bank)
        []
    """
    # 处理空值或无反应的情况
    if not rxn_string or rxn_string == '-':
        return []
    else:
        # 按分隔符分割反应ID字符串
        rxn_ids = [rxn_id.strip() for rxn_id in rxn_string.split(spliter) if rxn_id.strip()]
        # 获取每个反应ID的详细信息，过滤掉None值
        RXN_details = [get_rxn_detail(rxn_id, rxn_bank) for rxn_id in rxn_ids]
        # 过滤掉None值，只返回有效的Reaction对象
        return [rxn for rxn in RXN_details if rxn is not None] 


def get_rxn_details_batch(df_rxns, rxn_bank, rxn_id_column='RXNRECer', spliter=cfg.SPLITER):
    """
    批量处理DataFrame中的反应数据，为每行添加反应详情信息
    
    Args:
        df_rxns (pd.DataFrame): 包含反应数据的DataFrame
        rxn_bank (pd.DataFrame): 反应数据库DataFrame
        rxn_id_column (str, optional): 反应ID列名，默认为'RXNRECer'
        spliter (str, optional): 分隔符，默认为配置文件中的SPLITER
        
    Returns:
        pd.DataFrame: 扩展后的DataFrame，新增列：
            - RXN_details: 包含该行所有Reaction对象的列表
            
    Examples:
        >>> df = pd.DataFrame({'RXNRECer': ['RHEA:12345', 'RHEA:67890|RHEA:11111']})
        >>> result = get_rxn_details_batch(df, rxn_bank)
        >>> result['RXN_details'][0]  # 第一行的反应详情
        [<Reaction object>]
        >>> result['RXN_details'][1]  # 第二行的反应详情（包含2个反应）
        [<Reaction object>, <Reaction object>]
        
    Note:
        此函数会为每行创建一个RXN_details列，包含该行所有Reaction对象的列表
        如果某行包含多个反应ID（用分隔符分隔），会返回包含多个Reaction对象的列表
        无效的反应ID会被过滤掉，不会出现在结果中
    """
    
    # Create a copy to avoid modifying original
    result_df = df_rxns.copy()
    pandarallel.initialize()
    result_df['RXN_details'] = result_df[rxn_id_column].parallel_apply(
        lambda x: get_rxn_details_list(x, rxn_bank, spliter)
    )
    
    return result_df


def merge_reaction_with_s3_info(RXN_details, s3_info):
    """
    将S3信息补充回每个反应中，生成便于前端解析的JSON数据
    
    Args:
        RXN_details (list): 反应详情列表，每个元素是Reaction对象
        s3_info (list): S3信息列表，每个元素包含reaction_id, selected, rank, confidence, reason
        
    Returns:
        list: 合并后的反应信息列表，每个元素是包含S3信息的字典
        
    Examples:
        >>> RXN_details = [reaction_obj1, reaction_obj2, reaction_obj3]
        >>> s3_info = [
        ...     {'reaction_id': 'RHEA:14709', 'selected': 'yes', 'rank': 1, 'confidence': 0.95, 'reason': '...'},
        ...     {'reaction_id': 'RHEA:24076', 'selected': 'no', 'confidence': 0.2, 'reason': '...'},
        ...     {'reaction_id': 'RHEA:32187', 'selected': 'no', 'confidence': 0.1, 'reason': '...'}
        ... ]
        >>> merged = merge_reaction_with_s3_info(RXN_details, s3_info)
        >>> # 结果包含完整的反应信息和S3评分信息
    """
    if not RXN_details or not s3_info:
        return []
    
    # 创建S3信息的查找字典，以reaction_id为key
    s3_lookup = {}
    for s3_item in s3_info:
        if 'reaction_id' in s3_item:
            s3_lookup[s3_item['reaction_id']] = s3_item
    
    merged_reactions = []
    
    for reaction in RXN_details:
        if reaction is None:
            continue
            
        # 获取反应的基本信息
        reaction_dict = reaction.to_dict()
        
        # 查找对应的S3信息
        s3_data = s3_lookup.get(reaction.reaction_id, {})
        
        # 合并信息
        enriched_reaction = {
            # 反应基本信息
            **reaction_dict,
            
            # S3评分信息
            's3_selected': s3_data.get('selected', 'no'),
            's3_rank': s3_data.get('rank', None),
            's3_confidence': s3_data.get('confidence', 0.0),
            's3_reason': s3_data.get('reason', ''),
            
            # 前端友好的字段
            'is_selected': s3_data.get('selected', 'no') == 'yes',
            'selection_rank': s3_data.get('rank', None),
            'confidence_score': s3_data.get('confidence', 0.0),
            'selection_reason': s3_data.get('reason', '')
        }
        
        merged_reactions.append(enriched_reaction)
    
    return merged_reactions


def create_frontend_friendly_json(RXN_details, s3_info, output_file=None):
    """
    创建前端友好的JSON文件，包含完整的反应信息和S3评分
    
    Args:
        RXN_details (list): 反应详情列表
        s3_info (list): S3信息列表
        output_file (str, optional): 输出文件路径，如果为None则返回字典
        
    Returns:
        dict or None: 如果output_file为None返回字典，否则返回None（写入文件）
        
    Examples:
        >>> # 生成JSON文件
        >>> create_frontend_friendly_json(RXN_details, s3_info, 'output.json')
        
        >>> # 返回字典数据
        >>> data = create_frontend_friendly_json(RXN_details, s3_info)
        >>> print(json.dumps(data, indent=2))
    """
    import time
    
    merged_data = merge_reaction_with_s3_info(RXN_details, s3_info)
    
    # 创建前端友好的数据结构
    frontend_data = {
        'reactions': merged_data,
        'summary': {
            'total_reactions': len(merged_data),
            'selected_reactions': len([r for r in merged_data if r.get('is_selected', False)]),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'data_source': 'RXNRECer S3 Analysis'
        }
    }
    
    if output_file:
        # 写入文件
        ftool.write_json_file(frontend_data, output_file)
        print(f"✅ 前端友好JSON文件已生成: {output_file}")
        return None
    else:
        # 返回字典数据
        return frontend_data
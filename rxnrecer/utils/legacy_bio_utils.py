'''
Author: Zhenkun Shi
Date: 2022-04-08 20:04:18
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-06-30 13:34:10
FilePath: /RXNRECer/tools/bioFunctionLib.py
Description: 序列比对方法

Copyright (c) 2022 by tibd, All Rights Reserved. 
'''

import pandas as pd
import os,sys,re, string, random
from datetime import datetime
import subprocess
from tqdm import tqdm
import requests
from requests.exceptions import RequestException
from tkinter import _flatten
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
from config import conf as cfg
from modules.structure.Tdi import Tdi
import tempfile
from Bio import SeqIO
from pathlib import Path
from itertools import chain
from tools import filetool 
import shlex
from itertools import combinations




# 读取 FASTA 文件并转换为 DataFrame
def fasta_to_dataframe(fasta_file):
    # 使用列表推导式简化数据提取
    data = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    
    # 创建 DataFrame
    df = pd.DataFrame(data, columns=["uniprot_id", "seq"])
    return df

#region 训练集、测试集获取序列比对结果（blast）
def getblast(train, test, k=1):
    """训练集、测试集获取序列比对结果（blast）
    Args:
        train (DataFrame): 参造数据  
        test (DataFrame): 需要比对的数据
        k (int): blast 获取的结果个数
    Returns:
        DataFrame: 比对结果
    """
    # 使用临时文件保存 fasta 格式数据
    with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_train, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.tsv') as res_blast, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.dmnd') as db_dmnd:
        
        # 将训练集和测试集转换为 fasta 格式
        table2fasta(train, fasta_train.name)
        table2fasta(test, fasta_test.name)

        # 构建命令
        cmd1 = ["diamond", "makedb", "--in", fasta_train.name, "-d", db_dmnd.name, "--quiet"]
        cmd2 = ["diamond", "blastp", "-d", db_dmnd.name, "-q", fasta_test.name, "-o", res_blast.name, "-b5", "-c1", "-k", str(k), "--quiet"]
        
        # 运行命令
        subprocess.run(cmd1, check=True)
        subprocess.run(cmd2, check=True)

        # 读取比对结果
        res_data = pd.read_csv(res_blast.name, sep='\t', names=['id', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    
    return res_data
#endregion


#region 分隔底物产物字符串为化合物列表
def split_compound_str(componds_str, remove_reaction_coefficient=False):
    """ 分隔底物产物字符串为化合物列表

    Args:
        componds_str (string): 反应物或产物字符串
        remove_reaction_coefficient (bool, optional): 是否去除反应系数. Defaults to False.

    Returns:
        list: 反应物或产物列表
    """
    # 分割化合物
    list_componds = componds_str.split('#')
    
    # 去除反应系数
    if remove_reaction_coefficient:
        pattern = re.compile('(^[\d]{0,} )|(^n{0,1} )|(^a{0,1} )|(^an{0,1} )|(^[\d]/[\d]{0,1} )') # 数字系数、n系数、a系数、an系数、分数系数
        list_componds = [pattern.sub('', item) for item in list_componds]

    return list_componds

#endregion

#region 读取cdhit聚类结果
def get_cdhit_results(cdhit_clstr_file):
    """读取cdhit聚类结果

    Args:
        cdhit_clstr_file (string): 聚类结果文件

    Returns:
        DataFrame: ['cluster_id','uniprot_id','identity'， 'is_representative']
    """
    counter = 0
    res = []
    with open(cdhit_clstr_file,'r') as f:
        for line in f:
            if 'Cluster' in line:
                cluster_id = line.replace('>Cluster','').replace('\n', '').strip()
                continue
            str_uids= line.replace('\n','').split('>')[1].replace('at ','').split('... ')
                        
            if '*' in str_uids[1]:
                identity = 1
                isrep = True
            else:
                identity = float(str_uids[1].strip('%')) /100
                isrep = False

            res = res +[[cluster_id, str_uids[0], identity, isrep ]]

    resdf = pd.DataFrame(res, columns=['cluster_id','uniprot_id','identity', 'is_representative']) #转换为DataFrame
    return resdf
#endregion

def pycdhit(uniportid_seq_df, identity=0.4, thred_num=4):
    """CD-HIT 序列聚类

    Args:
        uniportid_seq_df (DataFrame): [uniprot_id, seq] 蛋白DataFrame
        identity (float, optional): 聚类阈值. Defaults to 0.4.
        thred_num (int, optional): 聚类线程数. Defaults to 4.

    Returns:
        聚类结果 DataFrame: [cluster_id, uniprot_id, identity, is_representative, cluster_size]
    """
    if identity < 0.4:
        raise ValueError('Identity is too low, CD-HIT does not support identity less than 0.4')

    # Determine word size based on identity threshold
    word_size = {0.7: 5, 0.6: 4, 0.5: 3, 0.4: 2}.get(identity, 5)

    # Generate timestamp and temporary file paths
    time_stamp_str = datetime.now().strftime("%Y-%m-%d_%H_%M_%S_") + ''.join(random.sample(string.ascii_letters + string.digits, 16))
    cd_hit_fasta = f'{cfg.TEMP_DIR}cdhit_test_{time_stamp_str}.fasta'
    cd_hit_results = f'{cfg.TEMP_DIR}cdhit_results_{time_stamp_str}'
    cd_hit_cluster_res_file = f'{cfg.TEMP_DIR}cdhit_results_{time_stamp_str}.clstr'

    try:
        # Write the input sequences to a FASTA file
        table2fasta(uniportid_seq_df, cd_hit_fasta)

        # Execute CD-HIT clustering
        cmd = [
            'cd-hit', 
            '-i', cd_hit_fasta,
            '-o', cd_hit_results,
            '-c', str(identity),
            '-n', str(word_size),
            '-T', str(thred_num),
            '-M', '0',
            '-g', '1',
            '-sc', '1',
            '-sf', '1'
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Process the clustering result
        cdhit_cluster = get_cdhit_results(cdhit_clstr_file=cd_hit_cluster_res_file)

        # Get cluster sizes and merge them with the cluster dataframe
        cluster_size = cdhit_cluster.cluster_id.value_counts().reset_index(name='cluster_size')
        cluster_size.rename(columns={'index': 'cluster_id'}, inplace=True)
        cdhit_cluster = pd.merge(cdhit_cluster, cluster_size, on='cluster_id', how='left')

    except subprocess.CalledProcessError as e:
        print(f"Error during CD-HIT execution: {e}")
        raise

    finally:
        # Clean up temporary files
        for temp_file in [cd_hit_fasta, cd_hit_results, cd_hit_cluster_res_file]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    return cdhit_cluster



#region 统计EC数
def stiatistic_ec_num(eclist):
    """统计EC数量

    Args:
        eclist (list): 可包含多功能酶的EC列表，用；分割

    Returns:
        int: 列表中包含的独立EC数量
    """
    eclist = list(eclist.flatten()) #展成1维
    eclist = _flatten([item.split(';') for item in eclist]) #分割多功能酶
    eclist = [item.strip() for item in eclist] # 去空格
    num_ecs = len(set(eclist))
    return num_ecs
#endregion


#region 创建foldseek数据库

def make_foldseek_db(prp_df, db_name):
    # Step 1: 定义目标目录
    dir_pdb = f'{cfg.DIR_FOLDSEEK_PDB}{db_name}/pdb/'
    dir_db = f'{cfg.DIR_FOLDSEEK_PDB}{db_name}/DB/{db_name}'

    # Step 2: 创建目录（如果不存在）
    if not os.path.exists(dir_pdb):
        os.makedirs(dir_pdb)
        
    if not os.path.exists(f'{cfg.DIR_FOLDSEEK_PDB}{db_name}/DB'):
        os.makedirs(f'{cfg.DIR_FOLDSEEK_PDB}{db_name}/DB')

    # Step 3: 复制 PDB 文件到目标目录
    print(f'Copying PDBs total:{len(prp_df)} files')
    prp_df.path_pdb.parallel_apply(lambda x: filetool.cp_pdb(src=x, dst=f'{dir_pdb}{os.path.basename(x)}'))

    # Step 4: 构建 FoldSeek 数据库创建命令
    foldseek_cmd = f'foldseek createdb {dir_pdb} {dir_db}'

    print(f"Executing command: {foldseek_cmd}")
    
    # Step 5: 执行命令并捕获输出
    result = subprocess.run(
        foldseek_cmd,
        shell=True,  # 必须启用 Shell 模式以解析字符串命令
        stdout=subprocess.PIPE,  # 捕获标准输出
        stderr=subprocess.PIPE,  # 捕获错误输出
        text=True  # 输出为文本
    )
    
    # # 打印输出以供调试
    # print("STDOUT:", result.stdout)
    # print("STDERR:", result.stderr)
    
    # 检查返回码
    if result.returncode != 0:
        raise RuntimeError(f"FoldSeek command failed with error: {result.stderr}")

    # 返回成功状态
    return 'success'
#endregion

#region 准备测试集PDB文件
def gather_test_pdb_db(prp_df, db_name):
    # Step 1: 定义目标目录
    
    dir_pdb_test = f'{cfg.DIR_FOLDSEEK_PDB}{db_name}/pdb_test/'
        
    if not os.path.exists(dir_pdb_test):
        os.makedirs(dir_pdb_test)

    # Step 3: 复制 PDB 文件到目标目录
    print(f'Copying PDBs total: {len(prp_df)} files, target directory: {dir_pdb_test}')
    prp_df.path_pdb.parallel_apply(lambda x: filetool.cp_pdb(src=x, dst=f'{dir_pdb_test}{os.path.basename(x)}'))

    print(f'Copying PDBs finished total:{len(prp_df)}')

    # 返回成功状态
    return 'success'
#endregion


#region 计算FoldSeek 3DI描述子
def get_fold_seek_3di(pdb_path, output_path=None, threads=40, verbosity=0):
    """
    Converts a PDB file to a 3DI descriptor using FoldSeek, removes the .dbtype file,
    and reads the resulting 3DI file into a Tdi object. Skips FoldSeek calculation if
    the 3DI file already exists.
    
    Parameters:
    - pdb_path (str): Path to the input PDB file.
    - output_path (str): Path to save the output 3DI file. If None, derives path based on cfg.DIR_DATASET_3DI and pdb_path.
    - threads (int): Number of threads to use for FoldSeek. Default is 10.
    - verbosity (int): Verbosity level for FoldSeek. Default is 0.
    
    Returns:
    - Tdi: Parsed 3DI descriptor as a Tdi object.
    """
    if not os.path.isfile(pdb_path):
        raise FileNotFoundError(f"The PDB file does not exist: {pdb_path}")

    # Default output path if not provided
    if output_path is None:
        output_path = os.path.join(
            cfg.DIR_DATASET_3DI,
            '/'.join(os.path.splitext(pdb_path)[0].split('/')[-2:]) + '.3di'
        )
    # print(f'{output_path}')
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Check if the 3DI file already exists
    if os.path.exists(output_path):
        if verbosity > 0:
            print(f"3DI file already exists. Skipping FoldSeek: {output_path}")
    else:
        # Construct the FoldSeek command
        command = [
            "foldseek", "structureto3didescriptor", 
            pdb_path, output_path, 
            "--threads", str(threads), 
            "-v", str(verbosity)
        ]
        # print(command)
        try:
            # Execute the command
            subprocess.run(command, capture_output=True, text=True, check=True)

            # Remove the .dbtype file if it exists
            dbtype_file = output_path + ".dbtype"
            if os.path.exists(dbtype_file):
                os.remove(dbtype_file)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"FoldSeek failed with error: {e.stderr.strip()}")

    # Read and return the 3DI file
    res_3di = Tdi()
    res_3di.read_3di_file(output_path)

    return res_3di
#endregion


def seqs2fasta(seq1_id, seq1, seq2_id, seq2, path1, path2):
    """将两条序列写入 fasta 文件"""
    with open(path1, 'w') as f1:
        f1.write(f'>{seq1_id}\n{seq1}\n')
    with open(path2, 'w') as f2:
        f2.write(f'>{seq2_id}\n{seq2}\n')

def blast2seq(seq1_id, seq1, seq2_id, seq2):
    """
    使用 BLASTP 对两条蛋白序列进行比对。
    
    Args:
        seq1_id (str): 第一个序列的 ID（作为数据库）
        seq1 (str): 第一个序列的氨基酸序列
        seq2_id (str): 第二个序列的 ID（作为查询）
        seq2 (str): 第二个序列的氨基酸序列

    Returns:
        pd.DataFrame: 比对结果（DataFrame 格式）
    """
    with tempfile.NamedTemporaryFile(suffix=".fasta") as fasta1, \
         tempfile.NamedTemporaryFile(suffix=".fasta") as fasta2, \
         tempfile.TemporaryDirectory() as tmpdir:

        # 写入 fasta
        seqs2fasta(seq1_id, seq1, seq2_id, seq2, fasta1.name, fasta2.name)

        # 构建数据库
        subprocess.run([
            "makeblastdb", "-in", fasta1.name, "-dbtype", "prot", "-out", f"{tmpdir}/db"
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # 运行 blastp
        output_file = f"{tmpdir}/blast.out"
        subprocess.run([
            "blastp", "-query", fasta2.name, "-db", f"{tmpdir}/db", "-outfmt", "6", "-out", output_file
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # 读取结果
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                   'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        df = pd.read_csv(output_file, sep="\t", names=columns)

    return df



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

def get_best_pdb(uniprot_id: str, save_path=''):
    """
    获取 UniProt ID 对应的最佳 PDB 结构。
    优先使用 PDBe 最佳结构，若找不到则回退使用 AlphaFold v4 模型。
    返回: {
        "pdb_id": str,
        "source": "rcsb" | "alphafold",
        "resolution": float | None,
        "method": str,
        "ligands": list
    } 或 None（全部失败）
    """
    try:
        # Step 1: 查询 PDBe API
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
    except (RequestException, ValueError) as e:
        print(f"[ERROR] Request failed for {uniprot_id}: {e}")
        data = {}

    if data and uniprot_id in data:
        try:
            df = json_to_dataframe(data)
            best_pdb, best_chain = select_best_pdb(df)
            print(f"[INFO] Best PDB ID for {uniprot_id}: {best_pdb}")
            
            # 获取 PDB 元数据
            try:
                pdb_api = f"https://data.rcsb.org/rest/v1/core/entry/{best_pdb}"
                pdb_info = requests.get(pdb_api, timeout=10).json()
                main_info, ligands = parse_pdb_info(pdb_info)
            except Exception as e:
                print(f"[WARNING] Failed to parse metadata for {best_pdb}: {e}")
                main_info, ligands = {}, []

            # 下载结构
            file_path = os.path.join(save_path, f"RCSB_{uniprot_id}_{best_pdb}.pdb")
            download_pdb(best_pdb, file_path)

            return {
                "pdb_id": best_pdb,
                "source": "rcsb",
                "resolution": main_info.get("resolution"),
                "method": main_info.get("method", "X-ray/EM"),
                "ligands": ligands
            }

        except Exception as e:
            print(f"[ERROR] Failed to process RCSB PDB for {uniprot_id}: {e}")

    # Step 2: 回退使用 AlphaFold v4 模型
    try:
        af_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        file_path = os.path.join(save_path, f"AF_{uniprot_id}_F1_model_v4.pdb")
        response = requests.get(af_url, timeout=10)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print(f"[INFO] Fallback AlphaFold model downloaded for {uniprot_id}")

            return {
                "pdb_id": f"AF-{uniprot_id}-F1-model_v4",
                "source": "alphafold",
                "resolution": None,
                "method": "AlphaFold v4",
                "ligands": []
            }
        else:
            print(f"[WARNING] AlphaFold model not available for {uniprot_id}")
    except Exception as e:
        print(f"[ERROR] Failed to download AlphaFold model for {uniprot_id}: {e}")

    return None
    

def run_tmalign(pdb1, pdb2):
    """运行 TM-align 并提取结构比对结果（包括 TM-score 和 RMSD）"""
    result = subprocess.run(["TMalign", str(pdb1), str(pdb2)],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL,
                            text=True)

    tm1, tm2, rmsd = None, None, None

    for line in result.stdout.splitlines():
        line = line.strip()
        if line.startswith("TM-score=") and "Chain_1" in line:
            tm1 = float(line.split('=')[1].split()[0])
        elif line.startswith("TM-score=") and "Chain_2" in line:
            tm2 = float(line.split('=')[1].split()[0])
        elif "RMSD=" in line and "Seq_ID" in line:
            try:
                rmsd = float(line.split("RMSD=")[1].split(",")[0])
            except:
                pass

    if tm1 is not None and tm2 is not None:
        return {
            'tm_score_chain1': round(tm1, 6),
            'tm_score_chain2': round(tm2, 6),
            'tm_score_avg': round((tm1 + tm2) / 2, 6),
            'tm_score_max': round(max(tm1, tm2), 6),
            'rmsd_tmalign': round(rmsd, 3) if rmsd is not None else None
        }

    return None


def run_rmsd(pdb1: str, pdb2: str):
    """使用 PyMOL align 计算两个结构的 CA RMSD"""
    pdb1_safe = shlex.quote(pdb1)
    pdb2_safe = shlex.quote(pdb2)

    cmd = f"""
    pymol -c -q -d '
    load {pdb1_safe}, obj1;
    load {pdb2_safe}, obj2;
    align obj1 and name CA, obj2 and name CA;
    quit'
    """

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            return None

        for line in result.stdout.splitlines():
            if "RMSD =" in line:
                try:
                    rmsd = float(line.split("=")[1].strip().split()[0])
                    return round(rmsd, 6)
                except ValueError:
                    pass

    except Exception:
        return None

    return None
    
def align_all_structures(pdb_dir, output_tsv="tmalign_results.tsv"):
    """遍历并比对该目录下所有结构，输出 TM-score 表格"""
    pdb_dir = Path(pdb_dir)
    pdb_files = sorted([
        str(p.resolve()) 
        for p in chain(pdb_dir.glob("*.pdb"), pdb_dir.glob("*.cif"))
    ])
    pdb_pairs = pd.DataFrame(combinations(pdb_files, 2), columns=['pdb1', 'pdb2'])

    # 运行 TM-align
    tmalign_results = pdb_pairs.parallel_apply(
        lambda row: run_tmalign(row['pdb1'], row['pdb2']), axis=1
    )

    # 拆解成多列
    tmalign_df = pd.json_normalize(tmalign_results)
    pdb_pairs = pd.concat([pdb_pairs, tmalign_df], axis=1)

    # 运行 PyMOL RMSD
    pdb_pairs['rmsd_pymol'] = pdb_pairs.parallel_apply(
        lambda row: run_rmsd(row['pdb1'], row['pdb2']), axis=1
    )

    # 结构名称标准化
    pdb_pairs['pdb1'] = pdb_pairs['pdb1'].apply(lambda x: Path(x).stem.replace('-F1-model_v4', ''))
    pdb_pairs['pdb2'] = pdb_pairs['pdb2'].apply(lambda x: Path(x).stem.replace('-F1-model_v4', ''))

    # 清洗无效行
    valid_cols = ['tm_score_avg', 'rmsd_pymol']
    pdb_pairs = pdb_pairs.dropna(subset=valid_cols)

    # 排序并保存
    pdb_pairs = pdb_pairs.sort_values(valid_cols, ascending=[False, True]).reset_index(drop=True)
    ordered_cols = ['pdb1', 'pdb2'] + [c for c in pdb_pairs.columns if c not in ('pdb1', 'pdb2')]
    pdb_pairs = pdb_pairs[ordered_cols]
    pdb_pairs.to_csv(output_tsv, sep='\t', index=False)

    return pdb_pairs
'''
Author: Zhenkun Shi
Date: 2022-04-08 20:04:18
LastEditors: Zhenkun Shi
LastEditTime: 2023-05-25 09:15:12
FilePath: /preaction/utils/bioinfo/bioFunctionLib.py
Description: 序列比对方法

Copyright (c) 2022 by tibd, All Rights Reserved. 
'''

import pandas as pd
import os,sys,re, string, random
from datetime import datetime
import subprocess
from tqdm import tqdm
from tkinter import _flatten
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
from config import conf as cfg
from modules.structure.Tdi import Tdi
import tempfile
from Bio import SeqIO
from tools import filetool 


#region DataFrame表格转fasta文件
def table2fasta(table, file_out):
    """DataFrame表格转fasta文件, 输入两列，【序列名称，序列】

    Args:
        table (DataFrame): 包含序列名称、序列的DataFame
        file_out (_type_): 输出fasta文件路径
    """
    file = open(file_out, 'w')
    for index, row in table.iterrows():
        file.write(f'>{row.iloc[0]}\n')
        file.write(f'{ row.iloc[1]}\n')
    file.close()
    # print('Write finished')
#endregion


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
        聚类结果 DataFrame: [cluster_id,uniprot_id,identity,is_representative,cluster_size]
    """
    if identity>=0.7:
        word_size = 5
    elif identity>=0.6:
        word_size = 4
    elif identity >=0.5:
        word_size = 3
    elif identity >=0.4:
        word_size =2
    else:
        word_size = 5

    # 定义输入输出文件名


    
    time_stamp_str = datetime.now().strftime("%Y-%m-%d_%H_%M_%S_")+''.join(random.sample(string.ascii_letters + string.digits, 16))
    cd_hit_fasta = f'{cfg.TEMP_DIR}cdhit_test_{time_stamp_str}.fasta'
    cd_hit_results = f'{cfg.TEMP_DIR}cdhit_results_{time_stamp_str}'
    cd_hit_cluster_res_file =f'{cfg.TEMP_DIR}cdhit_results_{time_stamp_str}.clstr'

    # 写聚类fasta文件
    table2fasta(uniportid_seq_df, cd_hit_fasta)

    # cd-hit聚类
    cmd = f'cd-hit -i {cd_hit_fasta} -o {cd_hit_results} -c {identity} -n {word_size} -T {thred_num} -M 0 -g 1 -sc 1 -sf 1 > /dev/null 2>&1'
    os.system(cmd)
    cdhit_cluster = get_cdhit_results(cdhit_clstr_file=cd_hit_cluster_res_file)

    cluster_size = cdhit_cluster.cluster_id.value_counts()
    cluster_size = pd.DataFrame({'cluster_id':cluster_size.index,'cluster_size':cluster_size.values})
    cdhit_cluster = cdhit_cluster.merge(cluster_size, on='cluster_id', how='left')
    
    cmd = f'rm -f {cd_hit_fasta} {cd_hit_results} {cd_hit_cluster_res_file}'
    os.system(cmd)

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

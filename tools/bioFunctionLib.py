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
import tempfile
from Bio import SeqIO


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


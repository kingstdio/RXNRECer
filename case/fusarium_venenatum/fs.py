'''
Author: Zhenkun Shi kingstdio@gmail.com
Date: 2025-05-28 10:06:25
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-05-29 17:35:29
FilePath: /RXNRECer/case/fusarium_venenatum/fs.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''


import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')
from config import conf as cfg
import pandas as pd
import numpy as np
import subprocess
import shutil
from itertools import chain
import gzip
from tools import filetool 
from tools import btools
from collections import defaultdict
from itertools import combinations
from tqdm import tqdm
import plotly.graph_objects as go
import plotly.express as px
import shlex
from pathlib import Path


#region 1. 获取Brenda EC 到 UniProt ID 的映射字典
def get_dict_brenda_ec2uniprotid():
    """
    构建 EC 号到 UniProt ID 的映射字典（来自 BRENDA 数据）。

    数据源: brenda_reaction_uniprot_dataset.feather 文件，包含字段 'uniprot_id' 和 'ec'。

    操作流程：
    - 加载数据并移除缺失的 EC 号记录；
    - 将每条记录的 EC 列按分号分割（支持多个 EC）；
    - 构建 dict[ec] = list of unique uniprot_id。

    返回:
        defaultdict(list): 映射字典，每个 EC 号对应一个 UniProt ID 列表
    """
    # 加载指定字段
    data_brenda = pd.read_feather(f'{cfg.DATA_ROOT}brenda/brenda_reaction_uniprot_dataset.feather')[['uniprot_id', 'ec']]

    # 移除缺失 EC 的记录
    data_brenda = data_brenda.dropna(subset=['ec']).reset_index(drop=True)

    # 初始化字典：EC -> UniProt ID 列表
    dict_ec_uniprotID = defaultdict(list)

    # 遍历每一行，构建映射
    for _, row in data_brenda.iterrows():
        uniprot_id = row['uniprot_id']
        ec_list = [ec.strip() for ec in row['ec'].split(';') if ec.strip()]
        for ec in ec_list:
            if uniprot_id not in dict_ec_uniprotID[ec]:
                dict_ec_uniprotID[ec].append(uniprot_id)

    return dict_ec_uniprotID
#endregion

def get_dict_expasy_ec2uniprotid():
    """
    构建 EC 号到 UniProt ID 的映射字典（来自 Expasy 数据）。

    数据源: ec_expasy.feather 文件，字段 'ec' 和 'ref_swissprot'（多个 UniProt ID 用分号分隔）

    操作流程：
    - 过滤 ref_swissprot 不为空的数据；
    - 拆分每个 EC 行中多个 UniProt ID；
    - 构建 dict[ec] = list of unique uniprot_id。

    返回:
        defaultdict(list): EC -> UniProt ID 列表的映射
    """
    data_expasy = pd.read_feather(f'{cfg.DATA_ROOT}expasy/ec_expasy.feather')

    # 只保留有 UniProt ID 的条目
    data_expasy = data_expasy[data_expasy.ref_swissprot != ''].reset_index(drop=True)

    # 处理 UniProt ID 字符串（取每项第一个 ID）
    data_expasy['ref_swissprot'] = data_expasy['ref_swissprot'].apply(
        lambda x: cfg.SPLITER.join([item.split(',')[0].strip() for item in x.split(';')])
    )

    # 构建 EC -> UniProt ID 映射
    dict_ec_uniprotID = defaultdict(list)

    for _, row in data_expasy.iterrows():
        ec = row['ec'].strip()
        uniprot_ids = [uid.strip() for uid in row['ref_swissprot'].split(cfg.SPLITER) if uid.strip()]
        for uid in uniprot_ids:
            if uid not in dict_ec_uniprotID[ec]:
                dict_ec_uniprotID[ec].append(uid)

    return dict_ec_uniprotID

#region 2.复制文件到目标目录
def copy_file(source_file_path, destination_dir):
    """
    复制文件到目标目录。
    
    - 如果源文件不存在或目标文件已存在，则跳过；
    - 如果目标目录不存在，则自动创建；
    - 如果源文件为 .gz 格式，会在复制过程中自动解压，输出为去除 .gz 后缀的文件。

    参数:
        source_file_path (str or Path): 源文件路径
        destination_dir (str or Path): 目标文件夹路径

    返回:
        Path: 成功复制（或解压）后的目标文件路径；如果跳过则返回 None
    """
    source = Path(source_file_path)
    if not source.exists():
        print(f"跳过：源文件不存在 -> {source}")
        return None

    target_dir = Path(destination_dir)
    target_dir.mkdir(parents=True, exist_ok=True)  # 自动创建目标目录

    # 如果是 .gz 文件，解压复制为去掉 .gz 后缀的目标文件
    if source.suffix == '.gz':
        target_path = target_dir / source.stem
        if target_path.exists():
            # print(f"跳过：目标文件已存在 -> {target_path}")
            return None

        with gzip.open(source, 'rb') as f_in, open(target_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        # print(f"已解压复制：{source} -> {target_path}")
        return target_path
    else:
        target_path = target_dir / source.name
        if target_path.exists():
            # print(f"跳过：目标文件已存在 -> {target_path}")
            return None

        shutil.copy(str(source), str(target_path))
        # print(f"已复制：{source} -> {target_path}")
        return target_path
#endregion    

#region 3.准备示例结构数据
def prepCase4EC(case_ec, res_fus_ven, dict_ec_uniprotID):
    """
    为指定 EC 号准备结构案例文件。

    步骤包括：
    1. 从 BRENDA Ground Truth 中获取对应 UniProt ID 的蛋白；
    2. 使用 btools 工具查找其最佳结构（PDB），并复制到案例目录；
    3. 从 res_fus_ven 中筛选预测为该 EC 的蛋白 ID，复制其结构文件到同一目录。

    参数:
        case_ec (str): EC 编号（如 "1.1.1.1"）
        res_fus_ven (pd.DataFrame): 包含 'clean' 和 'input_id' 的 DataFrame，表示预测结果
        dict_ec_uniprotID (dict): EC 到 UniProt ID 的映射字典

    返回:
        None
    """
    # 获取 BRENDA 中该 EC 对应的 UniProt ID 列表
    case_ref_protein_brenda = dict_ec_uniprotID.get(case_ec, [])
    if not case_ref_protein_brenda:
        print(f"[警告] EC {case_ec} 在 BRENDA 中无对应 UniProt ID")
    
    # 准备案例工作目录路径（以 EC 号为文件夹名，"." 替换为 "_"）
    case_work_dir = f'{cfg.CASE_DIR}fusarium_venenatum/ec{case_ec.replace(".", "_")}/'

    # 复制 BRENDA 参考蛋白的结构文件（通过 btools 获取最佳 PDB）
    for uniprot_id in case_ref_protein_brenda:
        pdb_path = btools.get_best_pdb(uniprot_id=uniprot_id)
        copy_file(pdb_path, case_work_dir)

    # 复制预测为该 EC 的 Fusarium 蛋白结构文件
    matched_ids = res_fus_ven[res_fus_ven.clean == case_ec].input_id.tolist()
    for input_id in matched_ids:
        pdb_path = f'{cfg.DATA_ROOT}structure/pdb/ncbi/{input_id[:9]}/{input_id}.pdb'
        copy_file(pdb_path, case_work_dir)
#endregion


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

def get_pdb_align_by_ec(case_ec,res_fus_ven,map_dict='brenda', output_tsv="tmalign_results.tsv", onlinemode=True):
    case_work_dir = f'{cfg.CASE_DIR}fusarium_venenatum/ec{case_ec.replace(".", "_")}/'
    if onlinemode:
        if map_dict=='brenda':
            dict_ec_uniprotID =get_dict_brenda_ec2uniprotid()
        elif map_dict=='expasy':
            dict_ec_uniprotID = get_dict_expasy_ec2uniprotid()
        else:
            print('map_dict should be "brenda" or "expasy"')
            return
        prepCase4EC(case_ec=case_ec, res_fus_ven=res_fus_ven, dict_ec_uniprotID=dict_ec_uniprotID)
    df_pdb_align=align_all_structures(case_work_dir, output_tsv=output_tsv)
    return df_pdb_align


def plot_tm_score_discrete_heatmap(df_pdb_align, value_col='tm_score_avg'):
    """
    使用 Plotly 绘制 TM-score 离散分级热图（以 0.1 为间隔分色），适用于 TM-score 分布。

    参数:
        df_pdb_align (pd.DataFrame): 含 'pdb1', 'pdb2' 和 TM-score 列的 DataFrame
        value_col (str): TM-score 字段名（默认 'tm_score_avg'）
    """
    assert value_col in df_pdb_align.columns, f"{value_col} not found in DataFrame."

    proteins = sorted(set(df_pdb_align['pdb1']) | set(df_pdb_align['pdb2']))
    matrix = pd.DataFrame(np.nan, index=proteins, columns=proteins)

    for _, row in df_pdb_align.iterrows():
        p1, p2 = row['pdb1'], row['pdb2']
        val = row[value_col]
        matrix.loc[p1, p2] = val
        matrix.loc[p2, p1] = val
    
    # 自己对自己设为1
    for p in proteins:
        matrix.loc[p, p] = 1.0

    z = matrix.values
    colorscale = 'RdBu'

    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=matrix.columns,
        y=matrix.index,
        zmin=0,
        zmax=1,
        colorscale=colorscale,
        colorbar=dict(
            tickmode='array',
            title=value_col
        )
    ))

    fig.update_layout(
        title=f"Discrete TM-score Heatmap ({value_col})",
        xaxis_nticks=len(proteins),
        yaxis_nticks=len(proteins),
        width=max(700, len(proteins) * 20),
        height=max(600, len(proteins) * 20)
    )
    fig.show()
    return fig
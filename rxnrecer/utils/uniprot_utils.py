'''
Author: Zhenkun Shi
Date: 2023-06-21 14:55:14
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-07-01 00:00:00
FilePath: rxnrecer/utils/uniprot_utils.py
Description: UniProt utility functions (snapshot parsing)

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''
from __future__ import annotations

import gzip
import os
import re
import time
from typing import List, Tuple
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Use new config module
from rxnrecer.config import config as cfg



# region Read data from gzip
def read_file_from_gzip(file_in_path: str, file_out_path: str, extract_type: str, save_file_type: str = 'tsv') -> None:
    """
    Parse data from UniProt swiss format gzip snapshot and export.

    Args:
        file_in_path: Input gzip file path
        file_out_path: 输出文件路径（.tsv 或 .feather）
        extract_type: 提取类型，可选 'with_ec' | 'without_ec' | 'full'
        save_file_type: 保存类型，'tsv' 或 'feather'
    """
    if save_file_type not in {"tsv", "feather"}:
        raise ValueError("save_file_type must be 'tsv' or 'feather'")

    # feather 模式先导出到中间 TSV，再转换
    if save_file_type == 'feather':
        outpath = file_out_path
        file_out_path = os.path.join(cfg.TEMP_DIR, 'temprecords.tsv')

    table_head = [
        'id', 'name', 'isenzyme', 'isMultiFunctional', 'functionCounts', 'ec_number',
        'ec_specific_level', 'date_integraged', 'date_sequence_update',
        'date_annotation_update', 'seq', 'seqlength'
    ]
    os.makedirs(os.path.dirname(file_out_path) or '.', exist_ok=True)
    try:
        counter = 0
        saver = 0
        with open(file_out_path, 'w') as file_write_obj:
            file_write_obj.writelines('\t'.join(table_head) + '\n')
            with gzip.open(file_in_path, 'rt') as handle:
                for record in tqdm(SeqIO.parse(handle, 'swiss'), position=1, leave=True):
                    res = process_record(record, extract_type=extract_type)
                    counter += 1
                    if counter % 10 == 0:
                        file_write_obj.flush()
                    if res:
                        saver += 1
                        file_write_obj.writelines('\t'.join(map(str, res)) + '\n')
        if save_file_type == 'feather':
            indata = pd.read_csv(file_out_path, sep='\t')
            indata.to_feather(outpath)
    except Exception as e:
        raise RuntimeError(f"Failed to read gzip and export data: {e}")
 # endregion

#region EC 解析：从描述中拆分 EC 号
def extract_ec_list(description: str) -> List[str]:
    """
    从描述文本中提取所有 EC 号，返回扁平化列表。

    规则：匹配所有 "EC=xxx" 片段，并将片段内可能存在的逗号分隔多个 EC 进一步展开。
    """
    if not description:
        return []
    ec_matches = re.findall(r"EC=([0-9\-,.;]+)", description)
    if not ec_matches:
        return []
    # 将每个匹配片段按逗号展开
    ec_list: List[str] = []
    for part in ec_matches:
        for item in part.split(','):
            item = item.strip()
            if item:
                ec_list.append(item)
    return ec_list
# endregion

#region EC 汇总：拼接字符串、判断多功能、统计数量
def summarize_ec_list(ec_list: List[str]) -> Tuple[str, bool, int]:
    """
    根据 EC 列表返回：标准化字符串、多功能标记、数量。
    """
    if not ec_list:
        return '-', False, 0
    ec_str = ','.join(ec_list)
    is_multi = len(ec_list) > 1
    count = len(ec_list)
    return ec_str, is_multi, count
# endregion

#region EC 特异性层级：计算最具体层级（4 为最具体）
def compute_ec_specific_level(ec_list: List[str]) -> int:
    """
    依据 EC 列表计算特异性层级：对每个 EC 号按连字符数量估算层级，取最具体值。
    单一 EC 例如 "1.2.3.4" → 4；含 "-" 越多越不具体。
    """
    if not ec_list:
        return 0
    best = 0
    for ec in ec_list:
        level = 4 - ec.count('-')
        if level > best:
            best = level
    return best
# endregion

#region 提取单条含有EC号的数据
def process_record(record, extract_type: str = 'with_ec') -> List:
    """
    提取单条 UniProt 记录的关键信息，并根据提取类型筛选。

    Args:
        record: Bio.SeqRecord 对象，表示一条 UniProt 记录。
        extract_type: 提取类型，可选 'with_ec'（仅含EC号）、'without_ec'（不含EC号）、'full'（全部）。

    Returns:
        List: 提取的字段列表，或根据筛选规则返回空列表。
    """
    try:
        # 1. 解析描述信息，判断是否为酶（含EC号）
        description = record.description or ''
        is_enzyme = 'EC=' in description

        # 2. 初始化相关变量
        is_multi_functional = False  # 是否多功能酶
        function_counts = 0          # EC号数量
        ec_specific_level = 0        # EC号特异性层级
        ec = '-'                     # EC号字符串

        # 3. 提取 EC 及相关属性（拆分 → 汇总 → 计算层级）
        if is_enzyme:
            ec_list = extract_ec_list(description)
            ec, is_multi_functional, function_counts = summarize_ec_list(ec_list)
            ec_specific_level = compute_ec_specific_level(ec_list)

        # 4. 提取其他字段
        id_ = str(record.id).strip()
        name = str(record.name).strip()
        seq = str(record.seq).strip()
        date_integrated = str(record.annotations.get('date', '')).strip()
        date_sequence_update = str(record.annotations.get('date_last_sequence_update', '')).strip()
        date_annotation_update = str(record.annotations.get('date_last_annotation_update', '')).strip()
        seqlength = len(seq)

        # 5. 组装结果
        res = [
            id_, name, is_enzyme, is_multi_functional, function_counts, ec,
            ec_specific_level, date_integrated, date_sequence_update,
            date_annotation_update, seq, seqlength
        ]

        # 6. 根据提取类型筛选返回
        if extract_type == 'full':
            return res
        if extract_type == 'with_ec':
            return res if is_enzyme else []
        if extract_type == 'without_ec':
            return [] if is_enzyme else res
        return []
    except Exception as e:
        # 捕获异常，返回空列表
        return []
# endregion

def run_exact_task(infile: str, outfile: str) -> None:
    start = time.process_time()
    extract_type = 'full'
    read_file_from_gzip(file_in_path=infile, file_out_path=outfile, extract_type=extract_type)
    end = time.process_time()
    print('finished use time %6.3f s' % (end - start))



if __name__ =="__main__":
    print('success')
    # start =  time.process_time()
    # in_filepath_sprot = cfg.FILE_LATEST_SPROT
    # out_filepath_sprot = cfg.FILE_LATEST_SPROT_FEATHER
    
    # in_filepath_trembl = cfg.FILE_LATEST_TREMBL
    # out_filepath_trembl = cfg.FILE_LATEST_TREMBL_FEATHER

    # extract_type ='full'


    # # read_file_from_gzip(file_in_path=in_filepath_sprot, file_out_path=out_filepath_sprot, extract_type=extract_type, save_file_type='feather')

    # read_file_from_gzip(file_in_path=in_filepath_trembl, file_out_path=out_filepath_trembl, extract_type=extract_type, save_file_type='feather')
    # end =  time.process_time()
    # print('finished use time %6.3f s' % (end - start))

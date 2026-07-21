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
        file_out_path: Output file path (.tsv or .feather)
        extract_type: Extraction type, options: 'with_ec' | 'without_ec' | 'full'
        save_file_type: Save type, 'tsv' or 'feather'
    """
    if save_file_type not in {"tsv", "feather"}:
        raise ValueError("save_file_type must be 'tsv' or 'feather'")


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


def extract_ec_list(description: str) -> List[str]:
    """Extract EC numbers from a UniProt description string."""
    if not description:
        return []
    ec_matches = re.findall(r"EC=([0-9\-,.;]+)", description)
    if not ec_matches:
        return []

    ec_list: List[str] = []
    for part in ec_matches:
        for item in part.split(','):
            item = item.strip()
            if item:
                ec_list.append(item)
    return ec_list
# endregion


def summarize_ec_list(ec_list: List[str]) -> Tuple[str, bool, int]:
    """Return a normalized EC string, multifunction flag, and EC count."""
    if not ec_list:
        return '-', False, 0
    ec_str = ','.join(ec_list)
    is_multi = len(ec_list) > 1
    count = len(ec_list)
    return ec_str, is_multi, count
# endregion


def compute_ec_specific_level(ec_list: List[str]) -> int:
    """Estimate the most specific EC hierarchy level in a list of EC numbers."""
    if not ec_list:
        return 0
    best = 0
    for ec in ec_list:
        level = 4 - ec.count('-')
        if level > best:
            best = level
    return best
# endregion


def process_record(record, extract_type: str = 'with_ec') -> List:
    """Extract selected fields from one UniProt SeqRecord."""
    try:

        description = record.description or ''
        is_enzyme = 'EC=' in description


        is_multi_functional = False
        function_counts = 0
        ec_specific_level = 0
        ec = '-'


        if is_enzyme:
            ec_list = extract_ec_list(description)
            ec, is_multi_functional, function_counts = summarize_ec_list(ec_list)
            ec_specific_level = compute_ec_specific_level(ec_list)


        id_ = str(record.id).strip()
        name = str(record.name).strip()
        seq = str(record.seq).strip()
        date_integrated = str(record.annotations.get('date', '')).strip()
        date_sequence_update = str(record.annotations.get('date_last_sequence_update', '')).strip()
        date_annotation_update = str(record.annotations.get('date_last_annotation_update', '')).strip()
        seqlength = len(seq)


        res = [
            id_, name, is_enzyme, is_multi_functional, function_counts, ec,
            ec_specific_level, date_integrated, date_sequence_update,
            date_annotation_update, seq, seqlength
        ]


        if extract_type == 'full':
            return res
        if extract_type == 'with_ec':
            return res if is_enzyme else []
        if extract_type == 'without_ec':
            return [] if is_enzyme else res
        return []
    except Exception as e:

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

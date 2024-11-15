# Standard Library Imports
import os
import sys

# Third-party Imports
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import plotly.graph_objects as go
from IPython.display import HTML
from pandarallel import pandarallel  # Importing pandarallel for parallel processing

# Setting up the path for the module
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1, '../')

# Local Imports
from config import conf as cfg
from tools import btools




# Read CSV files serially
def read_csv_files(file_paths):
    return [pd.read_csv(file, sep='\t') for file in file_paths]

# Function to get ec_rxn_nores
def get_ec_rxn_nores(pred_detail,  rxnkey):

    no_prediction = len(pred_detail[(pred_detail[rxnkey].str.contains('NO-PREDICTION'))])
    return [len(pred_detail), no_prediction]

def process_no_res(res_list,  rxnkey):
    return pd.DataFrame([get_ec_rxn_nores(pred_detail=res_list[item], rxnkey=rxnkey) for item in range(10)], 
                        columns=['test_size', 'no_prediction'])

# Make one-hot encoding label for each prediction
def make_labels(resdf, src_col1, src_col2, lb1, lb2, rxn_label_dict):
    resdf[[lb1, lb2]] = resdf.apply(
        lambda row: pd.Series({
            lb1: btools.make_label(reaction_id=str(row[src_col1]), rxn_label_dict=rxn_label_dict),
            lb2: btools.make_label(reaction_id=str(row[src_col2]), rxn_label_dict=rxn_label_dict)
        }), axis=1
    )
    return resdf

def apply_labels(res_list, src_col1, src_col2, lb1, lb2, rxn_label_dict):
    for i in tqdm(range(10)):
        res_list[i] = make_labels(resdf=res_list[i], src_col1=src_col1, src_col2=src_col2, lb1=lb1, lb2=lb2, rxn_label_dict=rxn_label_dict)
    return res_list


# Function to calculate metrics
def calculate_metrics(eva_df, ground_truth_col, pred_col, eva_name):
    res =  btools.rxn_eva_metric_with_colName(eva_df=eva_df, col_groundtruth=ground_truth_col, col_pred=pred_col, eva_name=eva_name)
    return res

# 多线程运行评价函数
def calculate_metrics_parallel(res_unirep, ground_truth_col, pred_col, max_workers=None):
    def run_metric_evaluation(index):
        return calculate_metrics(eva_df=res_unirep[index], ground_truth_col=ground_truth_col, pred_col=pred_col, eva_name=f'fold{index + 1}')
    
    results = [None] * len(res_unirep)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(run_metric_evaluation, i): i
            for i in range(len(res_unirep))
        }
        for future in as_completed(futures):
            i = futures[future]
            results[i] = future.result()
            
    results = pd.concat(results,axis=0).reset_index(drop=True)
    
    return results

def get_simi_Pred(pred_list, uniprot_rxn_dict, topk=3):
    uniprot_id_list = [item[0] for item in pred_list][:topk]
    rxn_ids = [uniprot_rxn_dict.get(uniprot_id) for uniprot_id in uniprot_id_list]
    rxn_res = (cfg.SPLITER).join(set(rxn_ids))
    return rxn_res

# Function to display results as HTML
def display_html_results(metrics, fold_std, eva_name):
    return HTML(f"""
         <div style="float:left; width:900px;">
              <h2 style='color:blue'>{eva_name} Evaluation 10 Fold Details</h2>
              {metrics.to_html()}
         </div>
         <div  style="float:left; width:600px;" >
              <h2 style='color:blue' >{eva_name} Evaluation 10 Fold Overview</h2>
                   {fold_std.to_html()}
         </div>
         """)
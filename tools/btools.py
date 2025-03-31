'''
Author: Zhenkun Shi
Date: 2023-04-20 06:23:40
LastEditors: Zhenkun Shi kingstdio@gmail.com
LastEditTime: 2025-03-24 18:53:50
FilePath: /preaction/pjlib/btools.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''

import json
import pandas as pd
import numpy as np
import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')
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
    for item in reaction_id.split(';'):
        array_index = rxn_label_dict.get(item)
        if array_index is not None:
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
    return ';'.join(ec_list)


if __name__ =='__main__':
    print('success')
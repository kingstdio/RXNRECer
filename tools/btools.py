'''
Author: Zhenkun Shi
Date: 2023-04-20 06:23:40
LastEditors: Zhenkun Shi
LastEditTime: 2023-05-17 09:44:00
FilePath: /preaction/pjlib/btools.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''

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
    
    
    # if reaction_id == 'EC-WITHOUT-REACTION':
    #     print('EC-WITHOUT-REACTION')
        
    # if reaction_id == 'NO-PREDICTION':
    #     print('NO-PREDICTION')
    
    resArray = [0]*len(rxn_label_dict)
    
    reactions = reaction_id.split(';')
    for item in reactions:
        try:
            array_index = rxn_label_dict.get(item)
            if array_index is not None:
                resArray[rxn_label_dict.get(item)] =1
            else:
                # resArray[1704]=1 # 单元素填错
                # print(reaction_id)
                resArray = [0]*len(rxn_label_dict)    #如果未找到平均分布到每个类中
        except Exception as e:
            print(reaction_id + '\n')
            
            
    
    return resArray
#endregion




#format EC resluts from clean
def load_clean_resluts(res_file):
    def format_ec(eclist):
        res = [item.split('/')[0].replace('EC:','')  for item in eclist if item!=None]
        res = ';'.join(res)
        return res
    
    with open(res_file, 'r') as f:
        data = [line.strip().split(',') for line in f.readlines()]
    data = pd.DataFrame(data)
    data['ec_clean'] = data.apply(lambda x : format_ec(eclist=x[1:].to_list()), axis=1)
    data=data[[0,'ec_clean']].rename(columns={0:'Entry'})
    return data

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
    return res_catfam


def load_ecpred_res(resfile):
    res_ecpred = pd.read_csv(f'{resfile}', sep='\t', header=0)
    res_ecpred = res_ecpred.rename(columns={'Protein ID':'id','EC Number':'ec_ecpred','Confidence Score(max 1.0)':'pident_ecpred'})
    res_ecpred['ec_ecpred']= res_ecpred.ec_ecpred.apply(lambda x : '-' if x=='non Enzyme' else x) 
    return res_ecpred

#endregion








if __name__ =='__main__':
    print('success')
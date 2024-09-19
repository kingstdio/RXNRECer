'''
Author: Zhenkun Shi
Date: 2023-03-16 11:11:37
LastEditors: Zhenkun Shi
LastEditTime: 2023-04-17 17:38:56
FilePath: /preaction/pjlib/commonfunction.py
Description: common function libs

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''

import re
import sys,os
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')

import pandas as pd
import numpy as np
import math
import pjlib.evaluation as eva

from tkinter import _flatten
from config import conf as cfg
#region ta
def ta_ds_spliter(dataset, test_ratio=0.2):
    
    test_size = math.ceil(len(dataset)*test_ratio) # 计算测试数据个数
    # 按3倍取测试数据
    if test_size*3 < len(dataset):
        sample_size = test_size*3
    else:
        sample_size = len(dataset)

    # 形成初步的训练集、测试集
    test_df = dataset.sample(n=sample_size) 
    train_df = dataset[~dataset.index.isin(test_df.index.values)]
    
    # 计算训练集中反应id，蛋白ID
    train_reaction_ids = set(train_df.reaction_id)
    train_uniprot_ids = set(train_df.uniprot_id)
    
    # 判断测试集中的反应和蛋白是否在训练集中出现过，没出现则剔出测试集
    pandarallel.initialize()
    test_df['ds_type']=test_df.parallel_apply(lambda x: 'test' if (x.reaction_id in train_reaction_ids) and (x.uniprot_id in train_uniprot_ids ) else 'train', axis=1)
    
    #二次采样，取足测试数据个数
    test_df = test_df.sample(n=test_size)

    #重新计算训练集与测试集
    train = dataset[~dataset.index.isin(test_df.index.values)]
    test = dataset[dataset.index.isin(test_df.index.values)]

    with pd.option_context('mode.chained_assignment', None):
        train['ds_type'] = 'train'
        test['ds_type'] = 'test'

    res_df = pd.concat([train, test], axis=0).sort_index()

    return res_df
#endregion

#region tb
def tb_ds_spliter(dataset, test_ratio=0.2):
    
    test_size = math.ceil(len(set(dataset.uniprot_id))*test_ratio) # 计算测试数据个数
    # 按3倍取测试数据
    if test_size*3 < len(dataset):
        sample_size = test_size*3
    else:
        sample_size = len(dataset)

    # 获取酶
    proteins = dataset['uniprot_id'].drop_duplicates().reset_index(drop=True)
    proteins = proteins.sample(n=sample_size)
    # 形成初步的训练集、测试集
    test_df = dataset[dataset.uniprot_id.isin(proteins.values)]
    train_df = dataset[~dataset.uniprot_id.isin(test_df.uniprot_id)]

    #训练集中的反应列表
    train_reaction_ids = list(set(train_df.reaction_id))
    # 判断测试集中的反应是否在训练集中出现过，没出现则剔出测试集
    with pd.option_context('mode.chained_assignment', None):
        pandarallel.initialize()
        test_df['ds_type']=test_df.parallel_apply(lambda x: 'test' if (x.reaction_id.strip() in train_reaction_ids) else 'train', axis=1)

    #二次采样，取足测试数据个数
    test_df = test_df[~test_df.uniprot_id.isin(test_df[test_df['ds_type']=='train'].uniprot_id)]
    proteins_2 = test_df['uniprot_id'].drop_duplicates().reset_index(drop=True)
    proteins_2 = proteins_2.sample(n=test_size)

    #重新计算训练集与测试集
    test = dataset[dataset.uniprot_id.isin(list(set(proteins_2.values)))]
    train = dataset[~dataset.uniprot_id.isin(list(set(proteins_2.values)))]

    with pd.option_context('mode.chained_assignment', None):
        train['ds_type'] = 'train'
        test['ds_type'] = 'test'

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df
#endregion

#region tc
def tc_ds_spliter(dataset, test_ratio=0.2):
    
    test_size = math.ceil(len(set(dataset.reaction_id))*test_ratio) # 计算测试数据个数
    # 按3倍取测试数据
    if test_size*3 < len(dataset):
        sample_size = test_size*3
    else:
        sample_size = len(dataset)

    # 获取反应
    reactions = dataset['reaction_id'].drop_duplicates().reset_index(drop=True)
    reactions = reactions.sample(n=sample_size)
    # 形成初步的训练集、测试集
    test_df = dataset[dataset.reaction_id.isin(reactions.values)]
    train_df = dataset[~dataset.reaction_id.isin(test_df.reaction_id)]

    #训练集中的酶列表
    train_protein_ids = list(set(train_df.uniprot_id))

    # 判断测试集中的反应是否在训练集中出现过，没出现则剔出测试集
    with pd.option_context('mode.chained_assignment', None):
        pandarallel.initialize()
        test_df['ds_type']=test_df.parallel_apply(lambda x: 'test' if (x.uniprot_id in train_protein_ids) else 'train', axis=1)

    #二次采样，取足测试数据个数
    test_df = test_df[~test_df.reaction_id.isin(test_df[test_df['ds_type']=='train'].reaction_id)]
    reactions_2 = test_df['reaction_id'].drop_duplicates().reset_index(drop=True)
    reactions_2 = reactions_2.sample(n=test_size)

    #重新计算训练集与测试集
    test = dataset[dataset.reaction_id.isin(list(set(reactions_2.values)))]
    train = dataset[~dataset.reaction_id.isin(list(set(reactions_2.values)))]

    with pd.option_context('mode.chained_assignment', None):
        train['ds_type'] = 'train'
        test['ds_type'] = 'test'

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df
#endregion

#region td
def td_ds_spliter(dataset, test_ratio=0.2, gama=0.85):
    
    test_size = math.ceil(len(set(dataset.reaction_id))*test_ratio) # 计算测试数据个数

    # 获取反应
    reactions = dataset['reaction_id'].drop_duplicates().reset_index(drop=True)
    reactions = reactions.sample(n=math.ceil(test_size*gama))
    
    # 形成初步的训练集、测试集
    test_df = dataset[dataset.reaction_id.isin(reactions.values)]
    train_df = dataset[~dataset.reaction_id.isin(test_df.reaction_id)]

    # 计算训练集中反应id，蛋白ID
    train_reaction_ids = set(train_df.reaction_id)
    train_uniprot_ids = set(train_df.uniprot_id)


    test = test_df[(~test_df.uniprot_id.isin(list(train_uniprot_ids)))&(~test_df.reaction_id.isin(list(train_reaction_ids)))]
    train = dataset[~dataset.index.isin(test.index)]


    with pd.option_context('mode.chained_assignment', None):
        train['ds_type'] = 'train'
        test['ds_type'] = 'test'

    res_df = pd.concat([train, test], axis=0).sort_index()
    return res_df
#endregion

#region stats for TDA-D
def show_stat_datasets(datasets, setname):

    train= datasets[datasets.ds_type=="train"]
    test= datasets[datasets.ds_type=="test"]

    protein_in_train = len(set(train.uniprot_id))
    reaction_in_train = len(set(train.reaction_id))

    protein_in_test = len(set(test.uniprot_id))
    reaction_in_test = len(set(test.reaction_id))


    print(60*'-')
    print(setname)
    print(60*'-')
    print(f'Trainning records:             {len(train)}')
    print(f'Proteins in trainning set:     {protein_in_train}')
    print(f'Reactions in trainning set:    {reaction_in_train}')
    print(f'Testing records:               {len(test)}')
    print(f'Proteins in testing set:       {protein_in_test}')
    print(f'Reactions in testing set:      {reaction_in_test}')
    print('')
#endregion






def format_train_valid_test(dataset, feature_reactions, feature_proteins, vali_ratio=0.2):

    dataset = dataset[['reaction_id','uniprot_id','label', 'ds_type']].reset_index(drop=True)
    dataset = dataset.merge(feature_reactions, on='reaction_id', how='left').merge(feature_proteins.rename(columns={'id':'uniprot_id'}), on='uniprot_id', how='left')
    
    train = dataset[dataset.ds_type=='train'].reset_index(drop=True)
    test = dataset[dataset.ds_type=='test'].reset_index(drop=True)

    vali = train.sample(frac=vali_ratio)
    train = train[~train.index.isin(vali.index)]

    vali= vali.reset_index(drop=True)
    train = train.reset_index(drop=True)

    train_X = train.iloc[:,4:]
    train_Y = train.label.to_list()

    vali_X = vali.iloc[:,4:]
    vali_Y = vali.label.to_list()

    test_X = test.iloc[:,4:]
    test_Y = test.label.to_list()

    return train_X, train_Y, vali_X, vali_Y, test_X, test_Y


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


    
# 评价酶非酶    
def eva_isenzyme(baselineName, res_df, category, print_flag=False):
    if category == 'ec':
        tp = res_df[(res_df[f'{category}_{baselineName}'].str.contains(r'\d', na=False)) & (res_df.isenzyme_groundtruth == True)].shape[0]
        tn = res_df[(res_df[f'{category}_{baselineName}'] == '-') & (res_df.isenzyme_groundtruth == False)].shape[0]
        fp = res_df[(res_df[f'{category}_{baselineName}'].str.contains(r'\d', na=False)) & (res_df.isenzyme_groundtruth == False)].shape[0]
        fn = res_df[(res_df[f'{category}_{baselineName}'] == '-') & (res_df.isenzyme_groundtruth == True)].shape[0]
    elif category == 'rxn':
        tp = res_df[(res_df[f'{category}_{baselineName}'].str.contains('RHEA', na=False)) & (res_df.rxn_groundtruth != '-')].shape[0]
        tn = res_df[(res_df[f'{category}_{baselineName}'] == '-') & (res_df.rxn_groundtruth == '-')].shape[0]
        fp = res_df[(res_df[f'{category}_{baselineName}'].str.contains('RHEA', na=False)) & (res_df.rxn_groundtruth == '-')].shape[0]
        fn = res_df[(res_df[f'{category}_{baselineName}'] == '-') & (res_df.rxn_groundtruth.str.contains('RHEA'))].shape[0]
    else:
        raise ValueError("Invalid category. Please choose 'ec' or 'rxn'.")

    up = res_df[(res_df[f'{category}_{baselineName}'] == 'NO-PREDICTION') & (res_df.isenzyme_groundtruth==True if category == 'ec' else res_df.rxn_groundtruth.str.contains('RHEA'))].shape[0]
    un = res_df[(res_df[f'{category}_{baselineName}'] == 'NO-PREDICTION') & (res_df.isenzyme_groundtruth==False if category == 'ec' else res_df.rxn_groundtruth == '-')].shape[0]

    acc = (tp + tn) / (tp + tn + fp + fn + up + un + 1.4E-45)
    precision = tp / (tp + fp + up + un + 1.4E-45)
    recall = tp / (tp + fn + up + un + 1.4E-45)
    f1 = 2 * precision * recall / (precision + recall + 1.4E-45)
    ppv = tp / (tp + fp + up + un + 1.4E-45)
    npv = tn / (tn + fn + un + un + 1.4E-45)

    if print_flag:
        print('{:<20} {:<14} {:<14} {:<15} {:<15} {:<20} {:<15} {:<10} {:<10} {:<10}{:<10}{:<10}{}'.format(
            'baselineName', 'accuracy', 'precision', 'recall', 'PPV(Sensitivity)', 'NPV(Specificity)', 'F1', 'TP', 'FP', 'FN', 'TN', 'UP', 'UN'))
        print('{:<20} {:<14.6f} {:<14.6f} {:<15.6f} {:<15.6f} {:<20.6f} {:<15.6f} {:<10} {:<10} {:<10}{:<10}{:<10}{}'.format(
            baselineName, acc, precision, recall, ppv, npv, f1, tp, fp, fn, tn, up, un))

    return [baselineName, acc, precision, recall, ppv, npv, f1, tp, fp, fn, tn, up, un]


if __name__ =='__main__':
    print('success')
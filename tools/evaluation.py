'''
Author: Zhenkun Shi
Date: 2023-03-17 04:54:02
LastEditors: Zhenkun Shi
LastEditTime: 2023-07-14 13:23:11
FilePath: /preaction/pjlib/evaluation.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''

import sys,os, math
sys.path.insert(0, os.path.dirname(os.path.realpath('__file__')))
sys.path.insert(1,'../')

from sklearn import metrics
from joblib import Parallel, delayed
import pandas as pd





def calculate_precision_at_k(y_true, y_pred, k):
    tp = 0
    for i, pred in enumerate(y_pred):
        if set(y_true[i]).issubset(set(pred[:k])):
            tp += 1
    precision = tp / (k * len(y_true))
    return precision

def calculate_ndcg_at_k(y_true, y_pred, k):
    dcg = 0.0
    idcg = 0.0
    for i, pred in enumerate(y_pred):
        if set(y_true[i]).issubset(set(pred[:k])):
            rel = 1.0
        else:
            rel = 0.0
        if i < k:
            dcg += rel / math.log2(i+2)
        idcg += 1.0 / math.log2(i+2)
    ndcg = dcg / idcg
    return ndcg

def calculate_recall_at_k(y_true, y_pred, k):
    tp = 0
    total_relevant = sum([len(true) for true in y_true])
    for i, pred in enumerate(y_pred):
        
        for j in range(k):
            if  (len(pred) > 0) and (j<len(pred)):
                if pred[j] in y_true[i]:
                    tp += 1
    recall = tp / total_relevant
    return recall

def calculate_top_k_accuracy(y_true, y_pred, k):
    
    correct = 0
    for i, pred in enumerate(y_pred):
        if set(y_true[i]).issubset(set(pred[:k])):
            correct += 1
    accuracy = correct / len(y_true)
    return accuracy

def calculate_map_at_k(y_true, y_pred, k):
    ap = 0.0
    for i, pred in enumerate(y_pred):
        tp = 0.0
        precisions = []
        for j in range(k):
            if  (len(pred) > 0) and (j<len(pred)):
                if pred[j] in y_true[i]:
                    tp += 1
                    precisions.append(tp / (j+1))
        if len(precisions) > 0:
            ap += sum(precisions) / len(precisions)
    map = ap / len(y_true)
    return map


def calculate_metrics_multi_joblib(groundtruth, predict, average_type, print_flag=False):
    results = Parallel(n_jobs=4, timeout=6000)(
        delayed(metric_fn)(groundtruth, predict, average_type) for metric_fn in [
            lambda gt, pr, avg: metrics.accuracy_score(gt, pr),
            lambda gt, pr, avg: metrics.precision_score(gt, pr, average=avg, zero_division=True),
            lambda gt, pr, avg: metrics.recall_score(gt, pr, average=avg, zero_division=True),
            lambda gt, pr, avg: metrics.f1_score(gt, pr, average=avg, zero_division=True)
        ]
    )

    acc, precision, recall, f1 = results

    if print_flag:
        print('%.6f ' % acc, '\t%.6f' % precision, '\t%.6f' % recall, '\t%.6f' % f1, '\t%12s' %average_type)
    
    return [acc, precision, recall, f1, average_type]



def caculateMetrix(groundtruth, predict, baselineName='', type='binary', averege_type='macro', print_flag=True):
    
    if type == 'binary':
        acc = metrics.accuracy_score(groundtruth, predict)
        precision = metrics.precision_score(groundtruth, predict, zero_division=True )
        recall = metrics.recall_score(groundtruth, predict,  zero_division=True)
        f1 = metrics.f1_score(groundtruth, predict, zero_division=True)
        tn, fp, fn, tp = metrics.confusion_matrix(groundtruth, predict).ravel()
        
        ppv = tp/(tp+fn+1.4E-45)
        npv = tn/(fn+tn+1.4E-45)

        # print('%12s'%baselineName, '\t\t%.6f' %acc,'\t%.6f'% precision,'\t%.6f'%npv,'\t%.6f'% recall,'\t%.6f'% f1, '\t', 'tp:',tp,'fp:',fp,'fn:',fn,'tn:',tn)
        # print('{:<24} {:<14.6f} {:<25.6f} {:<15.6f} {:<15.6f} {:<20.6f} tp: {} fp: {} fn: {} tn: {}'.format(baselineName, acc, precision, npv, recall, f1, tp, fp, fn, tn))
        if print_flag:
            print('{:<24} {:<14} {:<25} {:<15} {:<15} {:<20} {:<15} {:<15} {:<15} {}'.format('baselineName', 
                                                                                             'accuracy', 
                                                                                             'precision', 
                                                                                             'recall', 
                                                                                             'PPV(Sensitivity)', 
                                                                                             'NPV(Specificity)', 
                                                                                             'F1', 'TP', 'FP', 'FN', 'TN'))
            
            print('{:<24} {:<14.6f} {:<25.6f} {:<15.6f} {:<15.6f} {:<20.6f} {:<15} {:<15} {:<15} {}'.format(
                                                                                            baselineName,
                                                                                            acc, 
                                                                                            precision, 
                                                                                            recall,
                                                                                            ppv,
                                                                                            npv, 
                                                                                            f1, 
                                                                                            tp, fp, fn, tn))


        return [baselineName, acc, precision, recall, ppv, npv,  f1, tp, fp, fn, tn]
    
    if type =='include_unfind':
        evadf = pd.DataFrame()
        evadf['g'] = groundtruth
        evadf['p'] = predict

        evadf_hot = evadf[~evadf.p.isnull()]
        evadf_cold = evadf[evadf.p.isnull()]

        tp = len(evadf_hot[(evadf_hot.g.astype('int')==1) & (evadf_hot.p.astype('int')==1)])
        fp = len(evadf_hot[(evadf_hot.g.astype('int')==0) & (evadf_hot.p.astype('int')==1)])        
        tn = len(evadf_hot[(evadf_hot.g.astype('int')==0) & (evadf_hot.p.astype('int')==0)])
        fn = len(evadf_hot[(evadf_hot.g.astype('int')==1) & (evadf_hot.p.astype('int')==0)])
        up = len(evadf_cold[evadf_cold.g==1])
        un = len(evadf_cold[evadf_cold.g==0])
        acc = (tp+tn)/(tp+fp+tn+fn+up+un)
        precision = tp/(tp+fp)
        npv = tn/(tn+fn)
        recall = tp/(tp+fn+up)
        f1=(2*precision*recall)/(precision+recall)
        print(  baselineName, 
                '\t%.6f' %acc,
                '\t%.6f'% precision,
                '\t%.6f'%npv,
                '\t%.6f'% recall,
                '\t%.6f'% f1, '\t', 
                'tp:',tp,'fp:',fp,'fn:',fn,'tn:',tn, 'up:',up, 'un:',un)

        return [baselineName, acc, precision, npv, recall, f1, tp, fp, fn, tn, up, un]
    
    if type == 'multi':
        # acc = metrics.accuracy_score(groundtruth, predict)
        # precision = metrics.precision_score(groundtruth, predict, average=averege_type, zero_division=True )
        # recall = metrics.recall_score(groundtruth, predict, average=averege_type, zero_division=True)
        # f1 = metrics.f1_score(groundtruth, predict, average=averege_type, zero_division=True)
        # if print_flag:
        #     print('%12s'%baselineName, ' \t%.6f '%acc,'\t%.6f'% precision, '\t%.6f'% recall,'\t%.6f'% f1)
        # return [baselineName, acc, precision, recall, f1]
        
        res = calculate_metrics_multi_joblib(groundtruth=groundtruth, predict=predict, average_type=averege_type, print_flag=print_flag)
        return res
        
        




# 按反应groundtruth 展开反应df
def expand_df_with_reaction_ground_truth(df):
    expanded_list=[]
    for index, row in df.iterrows():
        # print(index, row.reaction_groundtruth)
        
        if len(row.reaction_groundtruth)<2:
            temprow = row.copy()
            temprow.reaction_groundtruth = row.reaction_groundtruth[0]
            expanded_list = expanded_list + [temprow.to_list()]
        else:
            for item in row.reaction_groundtruth:
                temprow = row.copy()
                temprow.reaction_groundtruth = item

                expanded_list = expanded_list + [temprow.to_list()]


    return pd.DataFrame(expanded_list, columns=df.columns)

def caculate_k_metrics_core(testdf, topk):
    accuracy = calculate_top_k_accuracy(y_true=testdf.reaction_groundtruth.to_list(), y_pred=testdf.reaction_pred, k=topk)
    precision = calculate_precision_at_k(y_true=testdf.reaction_groundtruth.to_list(), y_pred=testdf.reaction_pred, k=topk)
    ndcg = calculate_ndcg_at_k(y_true=testdf.reaction_groundtruth.to_list(), y_pred=testdf.reaction_pred, k=topk)
    recall = calculate_recall_at_k(y_true=testdf.reaction_groundtruth.to_list(), y_pred=testdf.reaction_pred, k=topk)
    map = calculate_map_at_k(y_true=testdf.reaction_groundtruth.to_list(), y_pred=testdf.reaction_pred, k=topk)

    return accuracy, precision, ndcg,recall, map

def caculate_k_metrics(test_groundtruth_predict):
    # test_groundtruth_predict = expand_df_with_reaction_ground_truth(df=test_groundtruth_predict)

    FLENGTH = 12
    res_head = f'{"":<{FLENGTH}} \t {"kAccuracy":<{FLENGTH}} \t  {"kPrecision":<{FLENGTH}} \t {"kNDCG":<{FLENGTH}} \t {"kRecall":<{FLENGTH}} \t {"kMAP":<{FLENGTH}}'
    print(res_head)

    res_str = ''
    TOPKS=[1,2,3,5,10,15,20,50, 100]
    for topk in TOPKS:
        accuracy, precision, ndcg,recall, map = caculate_k_metrics_core(testdf=test_groundtruth_predict, topk=topk)
        toplb = f'@top-k{topk}'
        # res_str = res_str + f'{toplb:<{11}}  {accuracy:12.6f} \t {precision:9.6f} \t {ndcg:8.6f} \t {recall:8.6f} \t {map:6.6f}\n'
        print(f'{toplb:<{11}}  {accuracy:12.6f} \t {precision:9.6f} \t {ndcg:8.6f} \t {recall:8.6f} \t {map:6.6f}')

    # return res_head + res_str


if __name__ == '__main__':
    print('success')
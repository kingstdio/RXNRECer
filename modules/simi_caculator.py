'''
Author: Zhenkun Shi
Date: 2023-07-28 14:37:09
LastEditors: Zhenkun Shi
LastEditTime: 2023-07-28 14:43:41
FilePath: /preaction/modules/simi_caculator.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''


import sys,os
sys.path.append(os.path.dirname(os.path.realpath('__file__')))
sys.path.append('../../')
import config.conf as cfg
import numpy as np
import pandas as pd
import concurrent.futures
from sklearn.metrics.pairwise import cosine_similarity # cos
from sklearn.metrics.pairwise import euclidean_distances # 欧氏距离
from pandarallel import pandarallel 

# 余弦相似性计算
def get_cosine_similarity(fx,fy):
    """余弦相似性

    Args:
        fx (_type_): _description_
        fy (_type_): _description_

    Returns:
        _type_: _description_
    """
    res = cosine_similarity(fx, fy)
    return res

# 欧式距离相似性计算
def get_euclidean_distances(fx, fy):
    """欧氏距离

    Args:
        fx (_type_): _description_
        fy (_type_): _description_

    Returns:
        _type_: _description_
    """
    res = euclidean_distances(fx, fy,Y_norm_squared=None, squared=False)
    return res

# # 皮尔逊相似性计算
# def get_pearsonr_smilarity(fx, fy):
#     res = pearsonr(fx, fy)
#     return res

if __name__ =='__main__':

    

    data = pd.read_feather(cfg.FILE_WEB_REACTIONS_FEATURE_RNXFP)
    pandarallel.initialize( progress_bar=True)
    simi_euclidean = data.parallel_apply(lambda x: np.round(get_euclidean_distances(fx=[x[1:].values], fy=data.iloc[:,1:].values)[0], 4), axis=1)
    simi_cosine = data.parallel_apply(lambda x: np.round(get_cosine_similarity(fx=[x[1:].values], fy=data.iloc[:,1:].values)[0], 4), axis=1)

 
    data['euclidean_simi'] = simi_euclidean
    data['cosine_simi'] = simi_cosine
    # data[['reaction_id', 'euclidean_simi', 'cosine_simi']].to_feather(cfg.FILE_WEB_REACTIONS_SIMI_RNXFP)


    # aa = get_pearsonr_smilarity(fx=x[0], fy=y[0])
    # print(aa)

    # pandarallel.initialize( progress_bar=True)
    
    # data['simi_eclidean']=data.equation_smiles.parallel_apply(lambda x: get_euclidean_distances(fx=[x], fy=data.equation_smiles.to_list()))

    # aa = get_euclidean_distances(fx=x, fy=y)
    # print(aa)
import numpy as np
from rxnrecer.config import config as cfg



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
    for item in reaction_id.split(cfg.SPLITER):
        array_index = rxn_label_dict.get(item)
        if array_index is not None:
            if not isinstance(array_index, int):
                print(f"[WARNING] Invalid index type: {type(array_index)} | value: {array_index} | item: {item}")
            resArray[array_index] = 1

    return resArray

#endregion

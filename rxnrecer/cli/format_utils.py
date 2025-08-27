import json
import numpy as np
from typing import Dict, List, Union, Any


def convert_numpy_types(obj):
    """
    递归转换numpy类型为Python原生类型，便于JSON序列化
    """
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj


def format_s3_output(RXNRECer_with_prob: Dict[str, float], 
                     RXNRECER_S3: List[Dict[str, Any]], 
                     RXN_details: List[Any]) -> Union[Dict[str, Any], List[Dict[str, Any]]]:
    """
    格式化S3输出，统一反应数据结构
    
    Args:
        RXNRECer_with_prob: 包含反应ID和置信度的字典
        RXNRECER_S3: S3预测结果列表
        RXN_details: 反应详情列表
    
    Returns:
        单个反应返回字典，多个反应返回字典列表
    """
    print("RXNRECer_with_prob:", RXNRECer_with_prob)
    print("RXNRECER_S3:", RXNRECER_S3)
    print("RXN_details:", RXN_details)
    
    # 检查是否只有一个元素且key为'-'（无催化活性）
    if len(RXNRECer_with_prob) == 1 and '-' in RXNRECer_with_prob:
        rxns = {
            'reaction_id': '-',
            'prediction_confidence': convert_numpy_types(RXNRECer_with_prob['-']),
            'reaction_details': {
                'reactants': [],
                'products': [],
                'reaction_ec': '-',
                'reaction_equation': '-',
                'reaction_equation_ref_chebi': '-',
                'reaction_smiles': '-',
            },
            'reaction_rxnrecer_s3': {
                'selected': RXNRECER_S3[0].get('selected', ''),
                'ranking': RXNRECER_S3[0].get('rank', 1),
                'confidence': convert_numpy_types(RXNRECER_S3[0].get('confidence', 0.0)),
                'reasoning': RXNRECER_S3[0].get('reason', ''),
            }
        }
        return rxns
    
    else:
        # 处理多个反应的情况
        rxns_list = []
        reaction_keys = list(RXNRECer_with_prob.keys())
        
        for i, item in enumerate(RXN_details):
            print(item.to_dict())
            # 获取对应的反应ID和置信度
            reaction_id = reaction_keys[i] if i < len(reaction_keys) else f"reaction_{i}"
            confidence = RXNRECer_with_prob.get(reaction_id, 0.0)
            # 转换为Python float并保留4位小数
            confidence = round(float(confidence), 4)
            
            # 获取对应的RXNRECER-S3信息
            s3_info = RXNRECER_S3[i] if i < len(RXNRECER_S3) else RXNRECER_S3[0]
            
            # 构建标准化的反应字典
            reaction_dict = {
                'reaction_id': reaction_id,
                'prediction_confidence': convert_numpy_types(confidence),
                'reaction_details': {
                    'reactants': [reactant.to_dict() for reactant in item.reactants],
                    'products': [product.to_dict() for product in item.products],
                    'reaction_ec': item.reaction_ec,
                    'reaction_equation': item.reaction_equation,
                    'reaction_equation_ref_chebi': item.reaction_equation_ref_chebi,
                    'reaction_smiles': item.reaction_smiles,
                },
                'reaction_rxnrecer_s3': {
                    'selected': s3_info.get('selected', ''),
                    'ranking': s3_info.get('rank', i + 1),
                    'confidence': convert_numpy_types(s3_info.get('confidence', 0.0)),
                    'reasoning': s3_info.get('reason', ''),
                }
            }
            rxns_list.append(reaction_dict)
        
        return rxns_list


def save_formatted_reactions(formatted_result: Union[Dict[str, Any], List[Dict[str, Any]]], 
                           filename: str = 'formatted_reactions.json') -> None:
    """
    保存格式化的反应结果到JSON文件
    
    Args:
        formatted_result: 格式化后的反应结果
        filename: 输出文件名
    """
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(formatted_result, f, indent=2, ensure_ascii=False)
    print(f"结果已保存到 {filename}")


def load_and_format_reactions(res_rxnrecer, sindex: int = 0) -> Union[Dict[str, Any], List[Dict[str, Any]]]:
    """
    从res_rxnrecer数据中加载并格式化反应结果
    
    Args:
        res_rxnrecer: RXNRECer结果数据
        sindex: 样本索引
    
    Returns:
        格式化后的反应结果
    """
    formatted_result = format_s3_output(
        res_rxnrecer.iloc[sindex].RXNRECer_with_prob, 
        res_rxnrecer.iloc[sindex]['RXNRECER-S3'], 
        res_rxnrecer.iloc[sindex].RXN_details
    )
    
    print("格式化结果:")
    print(json.dumps(formatted_result, indent=2, ensure_ascii=False))
    
    return formatted_result

from typing import Dict, List, Union, Any


def format_rxn_output(RXNRECer_with_prob: Dict[str, float], 
                          RXN_details: List[Any],
                          RXNRECER_S3: List[Dict[str, Any]] = None,
                          mode: str = 's3') -> List[Dict[str, Any]]:
    """
    格式化反应输出，统一反应数据结构
    
    Args:
        RXNRECer_with_prob: 包含反应ID和置信度的字典
        RXN_details: 反应详情列表
        RXNRECER_S3: S3预测结果列表（可选，仅S3模式需要）
        mode: 模式选择，'s2' 或 's3'
    
    Returns:
        反应字典列表
    """
    # 检查是否只有一个元素且key为'-'（无催化活性）
    if len(RXNRECer_with_prob) == 1 and '-' in RXNRECer_with_prob:
        rxns = {
            'reaction_id': '-',
            'prediction_confidence': round(RXNRECer_with_prob['-'], 4),
            'reaction_details': {
                'reactants': [],
                'products': [],
                'reaction_ec': '-',
                'reaction_equation': '-',
                'reaction_equation_ref_chebi': '-',
                'reaction_smiles': '-',
            }
        }
        
        # 如果是S3模式，添加额外的S3信息
        if mode == 's3' and RXNRECER_S3:
            rxns['reaction_rxnrecer_s3'] = {
                'selected': RXNRECER_S3[0].get('selected', ''),
                'ranking': RXNRECER_S3[0].get('rank', 1),
                'confidence': RXNRECER_S3[0].get('confidence', 0.0),
                'reasoning': RXNRECER_S3[0].get('reason', ''),
            }
        
        return [rxns]
    
    else:
        # 处理多个反应的情况
        rxns_list = []
        reaction_keys = list(RXNRECer_with_prob.keys())
        
        for i, item in enumerate(RXN_details):
            # 获取对应的反应ID和置信度
            reaction_id = reaction_keys[i] if i < len(reaction_keys) else f"reaction_{i}"
            confidence = round(float(RXNRECer_with_prob.get(reaction_id, 0.0)), 4)
            
            # 构建标准化的反应字典
            reaction_dict = {
                'reaction_id': reaction_id,
                'prediction_confidence': confidence,
                'reaction_details': {
                    'reactants': [reactant.to_dict() for reactant in item.reactants],
                    'products': [product.to_dict() for product in item.products],
                    'reaction_ec': item.reaction_ec,
                    'reaction_equation': item.reaction_equation,
                    'reaction_equation_ref_chebi': item.reaction_equation_ref_chebi,
                    'reaction_smiles': item.reaction_smiles,
                }
            }
            
            # 如果是S3模式，添加额外的S3信息
            if mode == 's3' and RXNRECER_S3:
                s3_info = RXNRECER_S3[i] if i < len(RXNRECER_S3) else RXNRECER_S3[0]
                reaction_dict['reaction_rxnrecer_s3'] = {
                    'selected': s3_info.get('selected', ''),
                    'ranking': s3_info.get('rank', i + 1),
                    'confidence': s3_info.get('confidence', 0.0),
                    'reasoning': s3_info.get('reason', ''),
                }
            
            rxns_list.append(reaction_dict)
        
        return rxns_list

import sys,os
sys.path.insert(0, f"{os.path.dirname(os.path.realpath('__file__'))}/../../../")
from rxnrecer.config import config as cfg
import pandas as pd
import Reaction as rxnTool
import json

def test_fixed_reaction():
    """测试修复后的Reaction类处理失败案例"""
    try:
        # 读取之前失败的数据
        if os.path.exists('middle_failed_reactions.json'):
            with open('middle_failed_reactions.json', 'r') as f:
                failed_data = json.load(f)
            
            print(f'测试之前失败的 {len(failed_data)} 个反应...')
            
            success_count = 0
            still_failed = []
            
            for i, failed in enumerate(failed_data):
                rxn_id = failed['reaction_id']
                try:
                    reaction = rxnTool.Reaction(
                        failed['smiles'], 
                        failed['equation'], 
                        failed['equation_chebi'], 
                        rxn_id=rxn_id, 
                        rxn_ec=failed['ec_number']
                    )
                    
                    print(f'成功修复: {rxn_id} - 反应物:{len(reaction.reactants)}, 产物:{len(reaction.products)}')
                    success_count += 1
                    
                    # 尝试保存JSON
                    json_path = f'{cfg.DIR_RXN_JSON}{rxn_id.replace(":", "_")}.json'
                    reaction.save_json_file(json_path)
                    
                except Exception as e:
                    print(f'仍然失败: {rxn_id} - 错误: {e}')
                    still_failed.append({
                        'reaction_id': rxn_id,
                        'error': str(e),
                        'original_error': failed['error']
                    })
            
            print(f'\n修复结果:')
            print(f'成功修复: {success_count} 个')
            print(f'仍然失败: {len(still_failed)} 个')
            
            if still_failed:
                print('\n仍然失败的反应:')
                for fail in still_failed:
                    print(f"- {fail['reaction_id']}: {fail['error']}")
            
            return success_count, still_failed
        else:
            print('没有找到失败记录文件')
            return 0, []
        
    except Exception as e:
        print(f'测试过程出错: {e}')
        return 0, []

if __name__ == "__main__":
    success_count, still_failed = test_fixed_reaction()

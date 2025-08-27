import sys,os
sys.path.insert(0, f"{os.path.dirname(os.path.realpath('__file__'))}/../../../")
from rxnrecer.config import config as cfg
import pandas as pd
import Reaction as rxnTool
import json
import traceback

def middle_test():
    """测试中间部分数据"""
    try:
        # 读取数据
        rxns = pd.read_feather(cfg.FILE_RHEA_REACTION)
        print(f'数据读取成功，共 {len(rxns)} 条记录')
        
        # 测试中间的数据：从5000到10000
        start_idx = 5000
        end_idx = 10000
        print(f'测试中间5000条记录 (索引 {start_idx} 到 {end_idx})...')
        
        failed_rxns = []
        success_count = 0
        
        for i in range(start_idx, end_idx):
            rxn_row = rxns.iloc[i]
            rxn_id = rxn_row['reaction_id']
            
            if (i - start_idx) % 500 == 0:
                print(f'处理进度: {i-start_idx}/{end_idx-start_idx} (成功: {success_count}, 失败: {len(failed_rxns)})')
            
            try:
                # 创建Reaction对象
                reaction = rxnTool.Reaction(
                    rxn_row['equation_smiles'], 
                    rxn_row['equation'], 
                    rxn_row['equation_chebi'], 
                    rxn_id=rxn_id, 
                    rxn_ec=rxn_row['ec_number']
                )
                
                # 保存JSON文件
                json_path = f'{cfg.DIR_RXN_JSON}{rxn_id.replace(":", "_")}.json'
                reaction.save_json_file(json_path)
                
                success_count += 1
                
            except Exception as e:
                error_msg = str(e)
                print(f'失败: {rxn_id} (索引{i}) - 错误: {error_msg}')
                failed_rxns.append({
                    'reaction_id': rxn_id,
                    'index': i,
                    'error': error_msg,
                    'error_type': type(e).__name__,
                    'smiles': rxn_row['equation_smiles'],
                    'equation': rxn_row['equation'],
                    'equation_chebi': rxn_row['equation_chebi'],
                    'ec_number': rxn_row['ec_number']
                })
                
                # 找到错误就保存并继续
                if len(failed_rxns) >= 20:
                    print(f'已找到{len(failed_rxns)}个错误，停止测试')
                    break
        
        # 输出结果
        print(f'\n测试完成!')
        print(f'成功: {success_count} 条')
        print(f'失败: {len(failed_rxns)} 条')
        
        if failed_rxns:
            print('\n失败的反应详情:')
            for failed in failed_rxns:
                print(f"- {failed['reaction_id']} (索引{failed['index']}): {failed['error']}")
            
            # 保存失败记录到文件
            failed_file = 'middle_failed_reactions.json'
            with open(failed_file, 'w') as f:
                json.dump(failed_rxns, f, indent=2, ensure_ascii=False)
            print(f'失败记录已保存到: {failed_file}')
        else:
            print('中间5000条记录都处理成功！')
        
        return failed_rxns
        
    except Exception as e:
        print(f'程序执行错误: {e}')
        traceback.print_exc()
        return []

if __name__ == "__main__":
    failed_rxns = middle_test()

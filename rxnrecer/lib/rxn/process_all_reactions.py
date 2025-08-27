import sys,os
sys.path.insert(0, f"{os.path.dirname(os.path.realpath('__file__'))}/../../../")
from rxnrecer.config import config as cfg
import pandas as pd
import Reaction as rxnTool
import json
import traceback
from datetime import datetime

def process_all_reactions():
    """处理所有反应数据并记录失败情况"""
    start_time = datetime.now()
    
    try:
        # 读取数据
        rxns = pd.read_feather(cfg.FILE_RHEA_REACTION)
        print(f'数据读取成功，共 {len(rxns)} 条记录')
        print(f'开始处理所有记录... ({start_time.strftime("%Y-%m-%d %H:%M:%S")})')
        
        failed_rxns = []
        success_count = 0
        warning_count = 0  # 有警告但成功的数量
        
        for i in range(len(rxns)):
            rxn_row = rxns.iloc[i]
            rxn_id = rxn_row['reaction_id']
            
            # 每1000条显示一次进度
            if i % 1000 == 0:
                elapsed = datetime.now() - start_time
                print(f'处理进度: {i}/{len(rxns)} (成功: {success_count}, 警告: {warning_count}, 失败: {len(failed_rxns)}) - 用时: {elapsed}')
            
            try:
                # 创建Reaction对象
                reaction = rxnTool.Reaction(
                    rxn_row['equation_smiles'], 
                    rxn_row['equation'], 
                    rxn_row['equation_chebi'], 
                    rxn_id=rxn_id, 
                    rxn_ec=rxn_row['ec_number']
                )
                
                # 检查是否有反应物和产物
                total_molecules = len(reaction.reactants) + len(reaction.products)
                if total_molecules == 0:
                    # 虽然创建成功但没有分子，记录为警告
                    warning_count += 1
                
                # 保存JSON文件
                json_path = f'{cfg.DIR_RXN_JSON}{rxn_id.replace(":", "_")}.json'
                reaction.save_json_file(json_path)
                
                success_count += 1
                
            except Exception as e:
                error_msg = str(e)
                print(f'失败: {rxn_id} (第{i+1}条) - 错误: {error_msg}')
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
        
        # 输出结果
        end_time = datetime.now()
        total_time = end_time - start_time
        
        print(f'\n处理完成! ({end_time.strftime("%Y-%m-%d %H:%M:%S")})')
        print(f'总用时: {total_time}')
        print(f'成功: {success_count} 条 ({success_count/len(rxns)*100:.2f}%)')
        print(f'有警告: {warning_count} 条 ({warning_count/len(rxns)*100:.2f}%)')
        print(f'失败: {len(failed_rxns)} 条 ({len(failed_rxns)/len(rxns)*100:.2f}%)')
        
        # 保存统计结果
        stats = {
            'start_time': start_time.isoformat(),
            'end_time': end_time.isoformat(),
            'total_time_seconds': total_time.total_seconds(),
            'total_records': len(rxns),
            'success_count': success_count,
            'warning_count': warning_count,
            'failed_count': len(failed_rxns),
            'success_rate': success_count/len(rxns)*100,
            'failed_rate': len(failed_rxns)/len(rxns)*100
        }
        
        with open('process_statistics.json', 'w') as f:
            json.dump(stats, f, indent=2, ensure_ascii=False)
        
        if failed_rxns:
            print('\n失败的反应详情:')
            for failed in failed_rxns[:10]:  # 只显示前10个
                print(f"- {failed['reaction_id']}: {failed['error']}")
            
            if len(failed_rxns) > 10:
                print(f"... 还有 {len(failed_rxns)-10} 个失败记录")
            
            # 保存失败记录到文件
            failed_file = 'all_failed_reactions.json'
            with open(failed_file, 'w') as f:
                json.dump(failed_rxns, f, indent=2, ensure_ascii=False)
            print(f'\n失败记录已保存到: {failed_file}')
            
            # 保存失败的反应ID列表
            failed_ids = [f['reaction_id'] for f in failed_rxns]
            with open('failed_reaction_ids.txt', 'w') as f:
                for rxn_id in failed_ids:
                    f.write(f'{rxn_id}\n')
            print(f'失败的反应ID列表已保存到: failed_reaction_ids.txt')
        else:
            print('所有记录都处理成功！')
        
        print(f'\n处理统计已保存到: process_statistics.json')
        return success_count, warning_count, failed_rxns
        
    except Exception as e:
        print(f'程序执行错误: {e}')
        traceback.print_exc()
        return 0, 0, []

if __name__ == "__main__":
    success_count, warning_count, failed_rxns = process_all_reactions()

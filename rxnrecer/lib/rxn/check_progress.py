import sys,os
sys.path.insert(0, f"{os.path.dirname(os.path.realpath('__file__'))}/../../../")
import json
import glob
from rxnrecer.config import config as cfg

def check_progress():
    """检查JSON文件生成进度"""
    
    # 检查生成的JSON文件数量
    json_pattern = f'{cfg.DIR_RXN_JSON}*.json'
    json_files = glob.glob(json_pattern)
    json_count = len(json_files)
    
    print(f'已生成的JSON文件数量: {json_count}')
    
    # 检查统计文件
    if os.path.exists('process_statistics.json'):
        with open('process_statistics.json', 'r') as f:
            stats = json.load(f)
        print(f'处理统计:')
        print(f'- 总记录数: {stats["total_records"]}')
        print(f'- 成功: {stats["success_count"]} ({stats["success_rate"]:.2f}%)')
        print(f'- 失败: {stats["failed_count"]} ({stats["failed_rate"]:.2f}%)')
        print(f'- 总用时: {stats["total_time_seconds"]:.2f} 秒')
    
    # 检查失败记录
    if os.path.exists('failed_reaction_ids.txt'):
        with open('failed_reaction_ids.txt', 'r') as f:
            failed_ids = f.read().strip().split('\n')
        print(f'\n失败的反应ID数量: {len(failed_ids)}')
        if failed_ids:
            print('失败的反应ID (前10个):')
            for i, rxn_id in enumerate(failed_ids[:10]):
                print(f'  {i+1}. {rxn_id}')
    
    return json_count

if __name__ == "__main__":
    check_progress()

#!/usr/bin/env python3
"""
RXNRECer 数据验证脚本
用于验证用户下载的数据和模型文件是否完整
"""

import os
import json
import hashlib
from pathlib import Path
import argparse

def calculate_md5(file_path):
    """计算文件的MD5哈希值"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def verify_file(file_path, expected_md5=None):
    """验证单个文件"""
    if not os.path.exists(file_path):
        return False, "File not found"
    
    file_size = os.path.getsize(file_path)
    actual_md5 = calculate_md5(file_path)
    
    if expected_md5 and actual_md5 != expected_md5:
        return False, f"MD5 mismatch: expected {expected_md5}, got {actual_md5}"
    
    return True, {
        "size": file_size,
        "md5": actual_md5,
        "size_mb": file_size / (1024*1024)
    }

def verify_data_structure(data_root):
    """验证数据目录结构"""
    required_files = [
        "rhea/rhea_reactions.feather",
        "dict/dict_id2rxn.json",
        "dict/dict_rxn2id.json",
        "dict/dict_rxnrecers3_prompt.json",
        "feature_bank/featureBank.feather",
        "sample/sample10.fasta"
    ]
    
    print("🔍 验证数据目录结构...")
    missing_files = []
    
    for file_path in required_files:
        full_path = os.path.join(data_root, file_path)
        if os.path.exists(full_path):
            size = os.path.getsize(full_path)
            print(f"   ✅ {file_path} ({size / (1024*1024):.1f} MB)")
        else:
            missing_files.append(file_path)
            print(f"   ❌ {file_path} - 缺失")
    
    return len(missing_files) == 0, missing_files

def verify_model_structure(ckpt_root):
    """验证模型目录结构"""
    required_files = [
        "rxnrecer/production_185846best.pth"
    ]
    
    print("🔍 验证模型目录结构...")
    missing_files = []
    
    for file_path in required_files:
        full_path = os.path.join(ckpt_root, file_path)
        if os.path.exists(full_path):
            size = os.path.getsize(full_path)
            print(f"   ✅ {file_path} ({size / (1024*1024):.1f} MB)")
        else:
            missing_files.append(file_path)
            print(f"   ❌ {file_path} - 缺失")
    
    return len(missing_files) == 0, missing_files

def main():
    parser = argparse.ArgumentParser(description="RXNRECer数据验证脚本")
    parser.add_argument("--data-root", default="~/.rxnrecer/data", help="数据目录路径")
    parser.add_argument("--ckpt-root", default="~/.rxnrecer/ckpt", help="模型目录路径")
    parser.add_argument("--manifest", help="数据清单文件路径")
    
    args = parser.parse_args()
    
    # 展开用户路径
    data_root = os.path.expanduser(args.data_root)
    ckpt_root = os.path.expanduser(args.ckpt_root)
    
    print("🚀 RXNRECer 数据验证工具")
    print("=" * 50)
    
    # 检查目录是否存在
    if not os.path.exists(data_root):
        print(f"❌ 数据目录不存在: {data_root}")
        print("💡 请先下载并解压数据文件")
        print("💡 或者指定正确的数据目录: --data-root /path/to/data")
        return 1
    
    if not os.path.exists(ckpt_root):
        print(f"❌ 模型目录不存在: {ckpt_root}")
        print("💡 请先下载并解压模型文件")
        print("💡 或者指定正确的模型目录: --ckpt-root /path/to/ckpt")
        return 1
    
    # 验证数据文件
    print(f"\n📁 数据目录: {data_root}")
    data_valid, data_missing = verify_data_structure(data_root)
    
    # 验证模型文件
    print(f"\n📁 模型目录: {ckpt_root}")
    model_valid, model_missing = verify_model_structure(ckpt_root)
    
    # 总结
    print("\n" + "=" * 50)
    print("📊 验证结果总结:")
    
    if data_valid:
        print("   ✅ 数据文件验证通过")
    else:
        print(f"   ❌ 数据文件验证失败，缺失 {len(data_missing)} 个文件")
        for file in data_missing:
            print(f"      - {file}")
    
    if model_valid:
        print("   ✅ 模型文件验证通过")
    else:
        print(f"   ❌ 模型文件验证失败，缺失 {len(model_missing)} 个文件")
        for file in model_missing:
            print(f"      - {file}")
    
    if data_valid and model_valid:
        print("\n🎉 所有文件验证通过！RXNRECer可以正常使用。")
        print("\n💡 下一步:")
        print("   1. 配置LLM API密钥（如果需要S3模式）")
        print("   2. 运行测试: rxnrecer -i ~/.rxnrecer/data/sample/sample10.fasta -o test.tsv -m s1")
        return 0
    else:
        print("\n❌ 文件验证失败，请检查下载和解压过程。")
        print("💡 参考 DATA_DOWNLOAD.md 重新下载文件。")
        return 1

if __name__ == "__main__":
    exit(main())

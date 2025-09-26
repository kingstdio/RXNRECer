#!/usr/bin/env python3
"""
RXNRECer data validation script
Used to verify if user downloaded data and model files are complete
"""

import os
import json
import hashlib
from pathlib import Path
import argparse

def calculate_md5(file_path):
    """Calculate MD5 hash of file"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def verify_file(file_path, expected_md5=None):
    """Verify single file"""
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
    """Verify data directory structure"""
    required_files = [
        "rhea/rhea_reactions.feather",
        "dict/dict_id2rxn.json",
        "dict/dict_rxn2id.json",
        "dict/dict_rxnrecers3_prompt.json",
        "feature_bank/featureBank.feather",
        "sample/sample10.fasta"
    ]
    
    print("🔍 Verify data directory structure...")
    missing_files = []
    
    for file_path in required_files:
        full_path = os.path.join(data_root, file_path)
        if os.path.exists(full_path):
            size = os.path.getsize(full_path)
            print(f"   ✅ {file_path} ({size / (1024*1024):.1f} MB)")
        else:
            missing_files.append(file_path)
            print(f"   ❌ {file_path} - Missing")
    
    return len(missing_files) == 0, missing_files

def verify_model_structure(ckpt_root):
    """Verify model directory structure"""
    required_files = [
        "rxnrecer/production_185846best.pth"
    ]
    
    print("🔍 Verifying model directory structure...")
    missing_files = []
    
    for file_path in required_files:
        full_path = os.path.join(ckpt_root, file_path)
        if os.path.exists(full_path):
            size = os.path.getsize(full_path)
            print(f"   ✅ {file_path} ({size / (1024*1024):.1f} MB)")
        else:
            missing_files.append(file_path)
            print(f"   ❌ {file_path} - Missing")
    
    return len(missing_files) == 0, missing_files

def main():
    parser = argparse.ArgumentParser(description="RXNRECer data verification script")
    parser.add_argument("--data-root", default="~/.rxnrecer/data", help="Data directory path")
    parser.add_argument("--ckpt-root", default="~/.rxnrecer/ckpt", help="Model directory path")
    parser.add_argument("--manifest", help="Data manifest file path")
    
    args = parser.parse_args()
    
    # Expand user paths
    data_root = os.path.expanduser(args.data_root)
    ckpt_root = os.path.expanduser(args.ckpt_root)
    
    print("🚀 RXNRECer Data Verification Tool")
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
    print("📊 Validation result总结:")
    
    if data_valid:
        print("   ✅ 数据文件Validation passed")
    else:
        print(f"   ❌ 数据文件Validation failed，缺失 {len(data_missing)} 个文件")
        for file in data_missing:
            print(f"      - {file}")
    
    if model_valid:
        print("   ✅ 模型文件Validation passed")
    else:
        print(f"   ❌ 模型文件Validation failed，缺失 {len(model_missing)} 个文件")
        for file in model_missing:
            print(f"      - {file}")
    
    if data_valid and model_valid:
        print("\n🎉 所有文件Validation passed！RXNRECer可以正常使用。")
        print("\n💡 下一步:")
        print("   1. 配置LLM API密钥（如果需要S3模式）")
        print("   2. 运行测试: rxnrecer -i ~/.rxnrecer/data/sample/sample10.fasta -o test.tsv -m s1")
        return 0
    else:
        print("\n❌ 文件Validation failed，请检查下载和解压过程。")
        print("💡 参考 DATA_DOWNLOAD.md 重新下载文件。")
        return 1

if __name__ == "__main__":
    exit(main())

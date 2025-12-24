#!/usr/bin/env python
"""
示例：如何使用数据库配置管理系统

运行前请先初始化数据库：
    python scripts/init_db_config.py
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from rxnrecer.config.db_config import DatabaseConfigManager
from rxnrecer.utils.db_utils import DatabaseConnection


def example_basic_usage():
    """示例1：基本使用"""
    print("=" * 60)
    print("示例1：基本配置查询")
    print("=" * 60)
    
    # 创建配置管理器（使用环境变量或默认SQLite）
    config_manager = DatabaseConfigManager()
    
    # 查询配置项
    data_root = config_manager.get_config('DATA_ROOT')
    cache_enabled = config_manager.get_config('CACHE_ENABLED', default=False)
    batch_size = config_manager.get_config('BATCH_SIZE', default=32)
    
    print(f"DATA_ROOT: {data_root}")
    print(f"CACHE_ENABLED: {cache_enabled}")
    print(f"BATCH_SIZE: {batch_size}")
    print()


def example_feature_configs():
    """示例2：按功能获取配置"""
    print("=" * 60)
    print("示例2：获取功能的配置项")
    print("=" * 60)
    
    config_manager = DatabaseConfigManager()
    
    # 获取预测功能的所有配置
    prediction_configs = config_manager.get_feature_configs('prediction')
    
    print("预测功能的配置项：")
    for key, value in prediction_configs.items():
        print(f"  {key}: {value}")
    print()


def example_all_features():
    """示例3：列出所有功能"""
    print("=" * 60)
    print("示例3：列出所有功能")
    print("=" * 60)
    
    config_manager = DatabaseConfigManager()
    
    # 获取所有功能
    all_features = config_manager.get_all_features()
    
    print("所有功能：")
    for feature in all_features:
        status = "✓ 启用" if feature.get('enabled') else "✗ 禁用"
        print(f"  {status} - {feature['feature_name']}: {feature.get('feature_description', '无描述')}")
    print()


def example_update_config():
    """示例4：更新配置"""
    print("=" * 60)
    print("示例4：更新配置项")
    print("=" * 60)
    
    config_manager = DatabaseConfigManager()
    
    # 获取当前值
    old_batch_size = config_manager.get_config('BATCH_SIZE', default=32)
    print(f"当前 BATCH_SIZE: {old_batch_size}")
    
    # 更新配置
    new_batch_size = 64
    config_manager.set_config(
        config_key='BATCH_SIZE',
        config_value=new_batch_size,
        config_type='int',
        description='Batch size for model inference',
        category='performance'
    )
    
    # 清除缓存以获取最新值
    config_manager.clear_cache()
    
    # 验证更新
    updated_batch_size = config_manager.get_config('BATCH_SIZE')
    print(f"更新后 BATCH_SIZE: {updated_batch_size}")
    
    # 恢复原值
    config_manager.set_config('BATCH_SIZE', old_batch_size, 'int', 'Batch size', 'performance')
    config_manager.clear_cache()
    print(f"已恢复为: {config_manager.get_config('BATCH_SIZE')}")
    print()


def example_custom_connection():
    """示例5：自定义数据库连接"""
    print("=" * 60)
    print("示例5：自定义数据库连接")
    print("=" * 60)
    
    # 使用SQLite（示例）
    db_conn = DatabaseConnection(
        db_type='sqlite',
        database='example_config.db'
    )
    
    config_manager = DatabaseConfigManager(db_connection=db_conn)
    
    # 确保表存在
    config_manager.ensure_tables()
    
    # 设置一个测试配置
    config_manager.set_config(
        'TEST_CONFIG',
        'test_value',
        'string',
        'Test configuration',
        'test'
    )
    
    # 查询配置
    test_value = config_manager.get_config('TEST_CONFIG')
    print(f"测试配置值: {test_value}")
    
    # 清理
    db_conn.disconnect()
    print("已断开数据库连接")
    print()


def main():
    """运行所有示例"""
    print("\n" + "=" * 60)
    print("RXNRECer 数据库配置管理示例")
    print("=" * 60 + "\n")
    
    try:
        example_basic_usage()
        example_feature_configs()
        example_all_features()
        example_update_config()
        example_custom_connection()
        
        print("=" * 60)
        print("所有示例运行完成！")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ 错误: {e}")
        import traceback
        traceback.print_exc()
        print("\n提示：请先运行 'python scripts/init_db_config.py' 初始化数据库")


if __name__ == '__main__':
    main()


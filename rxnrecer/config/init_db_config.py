"""
Initialize database configuration with default settings and feature lists
"""
import os
import sys
from pathlib import Path
from rxnrecer.utils.db_utils import DatabaseConnection, get_db_connection_from_env
from rxnrecer.config.db_config import DatabaseConfigManager
from rxnrecer.utils import file_utils


def init_default_configs(config_manager: DatabaseConfigManager):
    """Initialize default configuration items"""
    
    # Project paths
    project_root = str(file_utils.get_project_root())
    
    default_configs = [
        # Data directories
        ('DATA_ROOT', f'{project_root}/data/', 'string', 'Root directory for data files', 'paths'),
        ('CACHE_DIR', f'{project_root}/results/cache/', 'string', 'Cache directory', 'paths'),
        ('TEMP_DIR', f'{project_root}/temp/', 'string', 'Temporary directory', 'paths'),
        ('RESULTS_DIR', f'{project_root}/results/', 'string', 'Results directory', 'paths'),
        
        # Model paths
        ('MODEL_DIR', f'{project_root}/ckpt/rxnrecer/', 'string', 'Model directory', 'paths'),
        ('PROSTT5_DIR', f'{project_root}/ckpt/prostt5', 'string', 'ProSTT5 model directory', 'paths'),
        
        # Feature files
        ('FILE_PRODUCTION_FEATURES', f'{project_root}/data/feature_bank/featureBank.feather', 'string', 'Production features file', 'files'),
        ('FILE_PRODUCTION_FEATURES_T5', f'{project_root}/files/features/featureBank_t5.feather', 'string', 'T5 features file', 'files'),
        
        # RHEA files
        ('FILE_RHEA_REACTION', f'{project_root}/data/rhea/rhea_reactions.feather', 'string', 'RHEA reactions file', 'files'),
        ('FILE_DS_RHEA_REACTIONS', f'{project_root}/data/rhea/ds_rhea_reactions.feather', 'string', 'RHEA dataset file', 'files'),
        
        # Dictionary files
        ('DICT_UNIPROT_RHEA', f'{project_root}/data/dict/dict_uniprot_rhea.json', 'string', 'UniProt-RHEA dictionary', 'files'),
        ('DICT_EC_RHEA', f'{project_root}/data/dict/dict_ec_rhea.json', 'string', 'EC-RHEA dictionary', 'files'),
        
        # LLM API configuration
        ('LLM_API_URL', os.environ.get('LLM_API_URL', ''), 'string', 'LLM API URL', 'api'),
        ('LLM_API_KEY', os.environ.get('LLM_API_KEY', ''), 'string', 'LLM API Key', 'api'),
        
        # Model configuration
        ('CACHE_ENABLED', 'True', 'bool', 'Enable caching', 'model'),
        ('CODE_MODE', 'PRODUCTION', 'string', 'Code mode', 'model'),
        ('SPLITER', ';', 'string', 'Default separator', 'model'),
        
        # Performance settings
        ('BATCH_SIZE', '32', 'int', 'Batch size for inference', 'performance'),
        ('NUM_WORKERS', '4', 'int', 'Number of worker threads', 'performance'),
        ('DEVICE', 'auto', 'string', 'Device to use (auto, cpu, cuda)', 'performance'),
    ]
    
    print("Initializing default configuration items...")
    for config_key, config_value, config_type, description, category in default_configs:
        config_manager.set_config(config_key, config_value, config_type, description, category)
        print(f"  ✓ Set {config_key} = {config_value}")
    
    print(f"\nInitialized {len(default_configs)} configuration items")


def init_default_features(config_manager: DatabaseConfigManager):
    """Initialize default feature list and associations"""
    
    features = [
        ('prediction', 'Reaction prediction functionality', True),
        ('embedding', 'Protein embedding generation', True),
        ('feature_extraction', 'Feature extraction from sequences', True),
        ('cache_management', 'Cache management and optimization', True),
        ('llm_reasoning', 'LLM-based reasoning for reaction prediction', False),
        ('visualization', 'Reaction and compound visualization', True),
        ('data_download', 'Data download and update functionality', True),
        ('evaluation', 'Model evaluation and metrics', True),
    ]
    
    print("\nInitializing default features...")
    for feature_name, description, enabled in features:
        config_manager.add_feature(feature_name, description, enabled)
        print(f"  ✓ Added feature: {feature_name} ({'enabled' if enabled else 'disabled'})")
    
    # Associate configurations with features
    print("\nAssociating configurations with features...")
    
    feature_configs = {
        'prediction': [
            'MODEL_DIR', 'FILE_PRODUCTION_FEATURES', 'FILE_PRODUCTION_FEATURES_T5',
            'BATCH_SIZE', 'DEVICE', 'CACHE_ENABLED'
        ],
        'embedding': [
            'PROSTT5_DIR', 'DEVICE', 'NUM_WORKERS'
        ],
        'feature_extraction': [
            'FILE_PRODUCTION_FEATURES', 'FILE_PRODUCTION_FEATURES_T5', 'NUM_WORKERS'
        ],
        'cache_management': [
            'CACHE_DIR', 'CACHE_ENABLED', 'TEMP_DIR'
        ],
        'llm_reasoning': [
            'LLM_API_URL', 'LLM_API_KEY', 'CACHE_ENABLED'
        ],
        'visualization': [
            'DATA_ROOT', 'RESULTS_DIR'
        ],
        'data_download': [
            'DATA_ROOT', 'TEMP_DIR'
        ],
        'evaluation': [
            'RESULTS_DIR', 'DATA_ROOT', 'BATCH_SIZE'
        ],
    }
    
    for feature_name, config_keys in feature_configs.items():
        for config_key in config_keys:
            try:
                config_manager.associate_config_with_feature(feature_name, config_key)
                print(f"  ✓ Associated {config_key} with {feature_name}")
            except Exception as e:
                print(f"  ⚠ Warning: Could not associate {config_key} with {feature_name}: {e}")
    
    print(f"\nInitialized {len(features)} features with configuration associations")


def main():
    """Main initialization function"""
    print("=" * 60)
    print("RXNRECer Database Configuration Initialization")
    print("=" * 60)
    
    try:
        # Get database connection
        print("\nConnecting to database...")
        db_conn = get_db_connection_from_env()
        print(f"  ✓ Connected to {db_conn.db_type} database")
        
        # Create config manager
        config_manager = DatabaseConfigManager(db_connection=db_conn)
        
        # Ensure tables exist
        print("\nCreating database tables...")
        config_manager.ensure_tables()
        print("  ✓ Tables created/verified")
        
        # Initialize default configurations
        init_default_configs(config_manager)
        
        # Initialize default features
        init_default_features(config_manager)
        
        print("\n" + "=" * 60)
        print("Initialization completed successfully!")
        print("=" * 60)
        print("\nYou can now use DatabaseConfigManager to:")
        print("  - Query configuration items: config_manager.get_config('CONFIG_KEY')")
        print("  - Get feature configs: config_manager.get_feature_configs('feature_name')")
        print("  - List all features: config_manager.get_all_features()")
        
    except Exception as e:
        print(f"\n❌ Error during initialization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    finally:
        if 'db_conn' in locals():
            db_conn.disconnect()


if __name__ == '__main__':
    main()


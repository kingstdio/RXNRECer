"""
Base settings for RXNRECer
"""

import os
from pathlib import Path
from typing import Dict, Any


class Settings:
    """Base settings class for RXNRECer."""
    
    def __init__(self):
        # Project root directory
        self.PROJECT_ROOT = Path(__file__).parent.parent.parent
        
        # Environment
        self.ENVIRONMENT = os.getenv('RXNRECER_ENV', 'development')
        self.DEBUG = self.ENVIRONMENT == 'development'
        
        # Data directories
        self.DATA_ROOT = self.PROJECT_ROOT / 'data'
        self.RAW_DATA_DIR = self.DATA_ROOT / 'raw'
        self.PROCESSED_DATA_DIR = self.DATA_ROOT / 'processed'
        self.MODELS_DIR = self.DATA_ROOT / 'models'
        self.SAMPLES_DIR = self.DATA_ROOT / 'samples'
        
        # Results directories
        self.RESULTS_ROOT = self.PROJECT_ROOT / 'results'
        self.PREDICTIONS_DIR = self.RESULTS_ROOT / 'predictions'
        self.EVALUATIONS_DIR = self.RESULTS_ROOT / 'evaluations'
        self.LOGS_DIR = self.RESULTS_ROOT / 'logs'
        
        # Temporary directory
        self.TEMP_DIR = self.PROJECT_ROOT / 'temp'
        
        # Create directories if they don't exist
        self._create_directories()
        
        # File separators
        self.SPLITER = ';'
        
        # Logging
        self.LOG_LEVEL = os.getenv('LOG_LEVEL', 'INFO')
        self.LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        
        # Device configuration
        self.DEVICE = os.getenv('DEVICE', 'auto')  # 'auto', 'cpu', 'cuda'
        
    def _create_directories(self):
        """Create necessary directories."""
        directories = [
            self.RAW_DATA_DIR,
            self.PROCESSED_DATA_DIR,
            self.MODELS_DIR,
            self.SAMPLES_DIR,
            self.PREDICTIONS_DIR,
            self.EVALUATIONS_DIR,
            self.LOGS_DIR,
            self.TEMP_DIR
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    
    def get_model_path(self, model_name: str) -> Path:
        """Get path for a specific model."""
        return self.MODELS_DIR / model_name
    
    def get_data_path(self, data_type: str, filename: str) -> Path:
        """Get path for a specific data file."""
        if data_type == 'raw':
            return self.RAW_DATA_DIR / filename
        elif data_type == 'processed':
            return self.PROCESSED_DATA_DIR / filename
        elif data_type == 'samples':
            return self.SAMPLES_DIR / filename
        else:
            raise ValueError(f"Unknown data type: {data_type}")
    
    def get_result_path(self, result_type: str, filename: str) -> Path:
        """Get path for a specific result file."""
        if result_type == 'predictions':
            return self.PREDICTIONS_DIR / filename
        elif result_type == 'evaluations':
            return self.EVALUATIONS_DIR / filename
        elif result_type == 'logs':
            return self.LOGS_DIR / filename
        else:
            raise ValueError(f"Unknown result type: {result_type}")


# Global settings instance
settings = Settings()

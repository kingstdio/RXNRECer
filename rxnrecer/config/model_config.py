"""
Model configuration for RXNRECer
"""

import torch
from dataclasses import dataclass
from typing import Optional, Dict, Any


@dataclass
class ModelConfig:
    """Configuration for RXNRECer model."""
    
    # Model architecture
    esm_out_dim: int = 1280
    gru_h_dim: int = 512
    att_dim: int = 32
    dropout_rate: float = 0.2
    freeze_esm_layers: int = 32
    output_dimensions: int = 10479
    
    # Training parameters
    batch_size: int = 32
    learning_rate: float = 1e-4
    weight_decay: float = 1e-5
    num_epochs: int = 100
    early_stopping_patience: int = 10
    
    # Prediction parameters
    prediction_batch_size: int = 100
    top_k: int = 5
    confidence_threshold: float = 0.5
    
    # Device configuration
    device: Optional[str] = None
    
    def __post_init__(self):
        """Set device after initialization."""
        if self.device is None:
            self.device = self._get_device()
    
    def _get_device(self) -> str:
        """Get the best available device."""
        if torch.cuda.is_available():
            return "cuda"
        else:
            return "cpu"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            'esm_out_dim': self.esm_out_dim,
            'gru_h_dim': self.gru_h_dim,
            'att_dim': self.att_dim,
            'dropout_rate': self.dropout_rate,
            'freeze_esm_layers': self.freeze_esm_layers,
            'output_dimensions': self.output_dimensions,
            'batch_size': self.batch_size,
            'learning_rate': self.learning_rate,
            'weight_decay': self.weight_decay,
            'num_epochs': self.num_epochs,
            'early_stopping_patience': self.early_stopping_patience,
            'prediction_batch_size': self.prediction_batch_size,
            'top_k': self.top_k,
            'confidence_threshold': self.confidence_threshold,
            'device': self.device
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'ModelConfig':
        """Create config from dictionary."""
        return cls(**config_dict)


@dataclass
class ESMConfig:
    """Configuration for ESM model."""
    
    model_name: str = "esm2_t33_650M_UR50D"
    repr_layers: list = None
    max_length: int = 1022
    truncation: bool = True
    
    def __post_init__(self):
        """Set default repr_layers."""
        if self.repr_layers is None:
            self.repr_layers = [33]


@dataclass
class EnsembleConfig:
    """Configuration for ensemble methods."""
    
    methods: list = None
    weights: list = None
    voting_method: str = 'soft'  # 'hard' or 'soft'
    
    def __post_init__(self):
        """Set default methods and weights."""
        if self.methods is None:
            self.methods = ['esm', 't5', 'unirep']
        if self.weights is None:
            self.weights = [1.0] * len(self.methods)


# Default configurations
DEFAULT_MODEL_CONFIG = ModelConfig()
DEFAULT_ESM_CONFIG = ESMConfig()
DEFAULT_ENSEMBLE_CONFIG = EnsembleConfig()


def get_model_config(config_name: str = 'default') -> ModelConfig:
    """Get model configuration by name."""
    configs = {
        'default': DEFAULT_MODEL_CONFIG,
        'small': ModelConfig(
            esm_out_dim=1280,
            gru_h_dim=256,
            att_dim=16,
            dropout_rate=0.1,
            output_dimensions=10479
        ),
        'large': ModelConfig(
            esm_out_dim=1280,
            gru_h_dim=1024,
            att_dim=64,
            dropout_rate=0.3,
            output_dimensions=10479
        )
    }
    
    return configs.get(config_name, DEFAULT_MODEL_CONFIG)

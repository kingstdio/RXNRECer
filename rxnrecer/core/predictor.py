"""
Core predictor for RXNRECer
"""

import torch
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
from pathlib import Path
import logging

from ..config.settings import settings
from ..config.model_config import ModelConfig, get_model_config
from ..utils.file_utils import fasta_to_dataframe, save_dataframe, read_json_file
from ..utils.bio_utils import validate_sequence, clean_sequence
from ..models.bgru_model import BGRUModel
from ..models.esm_embedding import ESMEmbedding


class RXNRECerPredictor:
    """Main predictor class for RXNRECer."""
    
    def __init__(self, model_config: Optional[ModelConfig] = None,
                 model_path: Optional[str] = None):
        """
        Initialize RXNRECer predictor.
        
        Args:
            model_config: Model configuration
            model_path: Path to pre-trained model
        """
        self.logger = logging.getLogger(__name__)
        
        # Set model configuration
        self.model_config = model_config or get_model_config()
        
        # Initialize model
        self.model = None
        self.esm_embedding = None
        self.device = torch.device(self.model_config.device)
        
        # Load model if path provided
        if model_path:
            self.load_model(model_path)
        
        # Reaction dictionary
        self.reaction_dict = None
        self._load_reaction_dict()
    
    def _load_reaction_dict(self):
        """Load reaction dictionary."""
        try:
            dict_path = settings.get_data_path('processed', 'reaction_dict.json')
            if dict_path.exists():
                self.reaction_dict = read_json_file(str(dict_path))
            else:
                self.logger.warning("Reaction dictionary not found")
        except Exception as e:
            self.logger.error(f"Failed to load reaction dictionary: {e}")
    
    def load_model(self, model_path: str):
        """
        Load pre-trained model.
        
        Args:
            model_path: Path to model file
        """
        try:
            # Initialize ESM embedding
            self.esm_embedding = ESMEmbedding(
                device=self.device,
                freeze_layers=self.model_config.freeze_esm_layers
            )
            
            # Initialize BGRU model
            self.model = BGRUModel(
                input_dimensions=self.model_config.esm_out_dim,
                gru_h_size=self.model_config.gru_h_dim,
                attention_size=self.model_config.att_dim,
                dropout=self.model_config.dropout_rate,
                output_dimensions=self.model_config.output_dimensions,
                device=self.device
            )
            
            # Load model weights
            checkpoint = torch.load(model_path, map_location=self.device)
            self.model.load_state_dict(checkpoint['model_state_dict'])
            self.model.eval()
            
            self.logger.info(f"Model loaded successfully from {model_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to load model: {e}")
            raise
    
    def predict(self, input_data: Union[str, pd.DataFrame],
                batch_size: Optional[int] = None,
                top_k: Optional[int] = None,
                include_equations: bool = False) -> pd.DataFrame:
        """
        Predict reactions for input sequences.
        
        Args:
            input_data: Input FASTA file path or DataFrame
            batch_size: Batch size for prediction
            top_k: Number of top predictions to return
            include_equations: Whether to include reaction equations
            
        Returns:
            DataFrame with predictions
        """
        # Set default parameters
        batch_size = batch_size or self.model_config.prediction_batch_size
        top_k = top_k or self.model_config.top_k
        
        # Load input data
        if isinstance(input_data, str):
            df = fasta_to_dataframe(input_data)
        else:
            df = input_data.copy()
        
        # Validate sequences
        df = self._validate_and_clean_sequences(df)
        
        if df.empty:
            raise ValueError("No valid sequences found in input data")
        
        # Perform prediction
        predictions = self._predict_batches(df, batch_size, top_k)
        
        # Add reaction equations if requested
        if include_equations:
            predictions = self._add_reaction_equations(predictions)
        
        return predictions
    
    def _validate_and_clean_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and clean protein sequences."""
        valid_sequences = []
        
        for idx, row in df.iterrows():
            sequence = row['seq']
            
            # Clean sequence
            cleaned_seq = clean_sequence(sequence)
            
            # Validate sequence
            if validate_sequence(cleaned_seq):
                valid_sequences.append({
                    'uniprot_id': row['uniprot_id'],
                    'seq': cleaned_seq
                })
            else:
                self.logger.warning(f"Invalid sequence for {row['uniprot_id']}")
        
        return pd.DataFrame(valid_sequences)
    
    def _predict_batches(self, df: pd.DataFrame, batch_size: int, top_k: int) -> pd.DataFrame:
        """Perform prediction in batches."""
        if self.model is None:
            raise ValueError("Model not loaded. Call load_model() first.")
        
        predictions = []
        
        # Process in batches
        for i in range(0, len(df), batch_size):
            batch_df = df.iloc[i:i+batch_size]
            
            # Get ESM embeddings
            embeddings = self.esm_embedding(batch_df['seq'].tolist())
            
            # Get model predictions
            with torch.no_grad():
                logits = self.model(embeddings)
                probabilities = torch.softmax(logits, dim=1)
            
            # Get top-k predictions
            top_probs, top_indices = torch.topk(probabilities, k=top_k, dim=1)
            
            # Convert to predictions
            for j, (_, row) in enumerate(batch_df.iterrows()):
                pred_row = {'uniprot_id': row['uniprot_id']}
                
                # Add top-k predictions
                pred_reactions = []
                pred_probabilities = []
                
                for k in range(top_k):
                    reaction_idx = top_indices[j, k].item()
                    probability = top_probs[j, k].item()
                    
                    if self.reaction_dict and str(reaction_idx) in self.reaction_dict:
                        reaction_id = self.reaction_dict[str(reaction_idx)]
                    else:
                        reaction_id = f"RXN_{reaction_idx:06d}"
                    
                    pred_reactions.append(reaction_id)
                    pred_probabilities.append(probability)
                
                pred_row['top_predictions'] = ';'.join(pred_reactions)
                pred_row['top_probabilities'] = ';'.join([f"{p:.4f}" for p in pred_probabilities])
                pred_row['confidence'] = pred_probabilities[0]
                
                predictions.append(pred_row)
        
        return pd.DataFrame(predictions)
    
    def _add_reaction_equations(self, predictions: pd.DataFrame) -> pd.DataFrame:
        """Add reaction equations to predictions."""
        # TODO: Implement reaction equation lookup
        # This would involve loading reaction database and matching reaction IDs
        return predictions
    
    def predict_ensemble(self, input_data: Union[str, pd.DataFrame],
                        ensemble_config: Optional[Dict] = None) -> pd.DataFrame:
        """
        Perform ensemble prediction using multiple models.
        
        Args:
            input_data: Input data
            ensemble_config: Ensemble configuration
            
        Returns:
            Ensemble predictions
        """
        # TODO: Implement ensemble prediction
        # This would involve loading multiple models and combining their predictions
        pass
    
    def save_predictions(self, predictions: pd.DataFrame, 
                        output_file: str,
                        output_format: str = 'tsv'):
        """
        Save predictions to file.
        
        Args:
            predictions: Predictions DataFrame
            output_file: Output file path
            output_format: Output format
        """
        save_dataframe(predictions, output_file, output_format)
        self.logger.info(f"Predictions saved to {output_file}")
    
    def evaluate_predictions(self, predictions: pd.DataFrame,
                           ground_truth: pd.DataFrame) -> Dict[str, float]:
        """
        Evaluate prediction performance.
        
        Args:
            predictions: Predictions DataFrame
            ground_truth: Ground truth DataFrame
            
        Returns:
            Dictionary with evaluation metrics
        """
        from ..utils.evaluation_utils import calculate_reaction_metrics
        
        metrics = calculate_reaction_metrics(ground_truth, predictions)
        return metrics

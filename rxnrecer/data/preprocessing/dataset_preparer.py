"""
Dataset preparation utilities for RXNRECer
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging

from ...utils.file_utils import save_dataframe, load_dataframe
from ...utils.bio_utils import validate_sequence, clean_sequence


class DatasetPreparer:
    """Dataset preparation and preprocessing utilities."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def prepare_task_dataset(self, 
                           protein_data: pd.DataFrame,
                           reaction_data: pd.DataFrame,
                           output_path: str,
                           min_sequence_length: int = 10,
                           max_sequence_length: int = 2000) -> pd.DataFrame:
        """
        Prepare dataset for reaction prediction task.
        
        Args:
            protein_data: DataFrame with protein sequences
            reaction_data: DataFrame with reaction information
            output_path: Path to save prepared dataset
            min_sequence_length: Minimum sequence length
            max_sequence_length: Maximum sequence length
            
        Returns:
            Prepared dataset DataFrame
        """
        self.logger.info("Preparing task dataset...")
        
        # Filter sequences by length
        protein_data = protein_data[
            (protein_data['seq'].str.len() >= min_sequence_length) &
            (protein_data['seq'].str.len() <= max_sequence_length)
        ].copy()
        
        # Clean and validate sequences
        protein_data['seq'] = protein_data['seq'].apply(clean_sequence)
        protein_data = protein_data[protein_data['seq'].apply(validate_sequence)]
        
        # Merge with reaction data
        dataset = pd.merge(protein_data, reaction_data, on='uniprot_id', how='inner')
        
        # Add dataset metadata
        dataset['sequence_length'] = dataset['seq'].str.len()
        dataset['is_valid'] = True
        
        # Save prepared dataset
        save_dataframe(dataset, output_path, 'feather')
        
        self.logger.info(f"Dataset prepared: {len(dataset)} samples")
        return dataset
    
    def prepare_3di_dataset(self,
                           protein_data: pd.DataFrame,
                           structure_data: pd.DataFrame,
                           output_path: str) -> pd.DataFrame:
        """
        Prepare dataset with 3D structure information.
        
        Args:
            protein_data: DataFrame with protein sequences
            structure_data: DataFrame with structure information
            output_path: Path to save prepared dataset
            
        Returns:
            Prepared dataset with 3D structure information
        """
        self.logger.info("Preparing 3DI dataset...")
        
        # Merge protein and structure data
        dataset = pd.merge(protein_data, structure_data, on='uniprot_id', how='inner')
        
        # Add structure metadata
        dataset['has_structure'] = dataset['pdb_id'].notna()
        dataset['structure_resolution'] = dataset.get('resolution', np.nan)
        
        # Save prepared dataset
        save_dataframe(dataset, output_path, 'feather')
        
        self.logger.info(f"3DI dataset prepared: {len(dataset)} samples")
        return dataset
    
    def split_dataset(self,
                     dataset: pd.DataFrame,
                     train_ratio: float = 0.8,
                     val_ratio: float = 0.1,
                     test_ratio: float = 0.1,
                     stratify_col: Optional[str] = None,
                     random_state: int = 42) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Split dataset into train/validation/test sets.
        
        Args:
            dataset: Input dataset
            train_ratio: Training set ratio
            val_ratio: Validation set ratio
            test_ratio: Test set ratio
            stratify_col: Column to stratify on
            random_state: Random seed
            
        Returns:
            Tuple of (train_df, val_df, test_df)
        """
        from sklearn.model_selection import train_test_split
        
        # Validate ratios
        assert abs(train_ratio + val_ratio + test_ratio - 1.0) < 1e-6, "Ratios must sum to 1.0"
        
        # First split: train vs rest
        train_df, rest_df = train_test_split(
            dataset, 
            test_size=(val_ratio + test_ratio),
            stratify=dataset[stratify_col] if stratify_col else None,
            random_state=random_state
        )
        
        # Second split: validation vs test
        val_size = val_ratio / (val_ratio + test_ratio)
        val_df, test_df = train_test_split(
            rest_df,
            test_size=(1 - val_size),
            stratify=rest_df[stratify_col] if stratify_col else None,
            random_state=random_state
        )
        
        # Add split labels
        train_df = train_df.copy()
        val_df = val_df.copy()
        test_df = test_df.copy()
        
        train_df['split'] = 'train'
        val_df['split'] = 'val'
        test_df['split'] = 'test'
        
        self.logger.info(f"Dataset split: train={len(train_df)}, val={len(val_df)}, test={len(test_df)}")
        
        return train_df, val_df, test_df
    
    def create_reaction_mapping(self,
                               dataset: pd.DataFrame,
                               output_path: str) -> Dict[str, int]:
        """
        Create reaction ID to index mapping.
        
        Args:
            dataset: Dataset with reaction information
            output_path: Path to save mapping
            
        Returns:
            Reaction ID to index mapping
        """
        unique_reactions = dataset['reaction_id'].unique()
        reaction_mapping = {rxn_id: idx for idx, rxn_id in enumerate(unique_reactions)}
        
        # Save mapping
        import json
        with open(output_path, 'w') as f:
            json.dump(reaction_mapping, f, indent=2)
        
        self.logger.info(f"Reaction mapping created: {len(reaction_mapping)} reactions")
        return reaction_mapping
    
    def create_protein_mapping(self,
                              dataset: pd.DataFrame,
                              output_path: str) -> Dict[str, int]:
        """
        Create protein ID to index mapping.
        
        Args:
            dataset: Dataset with protein information
            output_path: Path to save mapping
            
        Returns:
            Protein ID to index mapping
        """
        unique_proteins = dataset['uniprot_id'].unique()
        protein_mapping = {prot_id: idx for idx, prot_id in enumerate(unique_proteins)}
        
        # Save mapping
        import json
        with open(output_path, 'w') as f:
            json.dump(protein_mapping, f, indent=2)
        
        self.logger.info(f"Protein mapping created: {len(protein_mapping)} proteins")
        return protein_mapping

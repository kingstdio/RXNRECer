"""
Embedding generation utilities for RXNRECer
"""

import pandas as pd
import numpy as np
import torch
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging
import tempfile
import os

from ...utils.file_utils import save_dataframe, load_dataframe


class EmbeddingGenerator:
    """Embedding generation utilities for different models."""
    
    def __init__(self, device: Optional[str] = None):
        self.logger = logging.getLogger(__name__)
        self.device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
    
    def generate_esm_embeddings(self,
                               dataset: pd.DataFrame,
                               output_path: str,
                               model_name: str = "esm2_t33_650M_UR50D",
                               layer: int = 33,
                               pooling: str = 'mean') -> pd.DataFrame:
        """
        Generate ESM embeddings for protein sequences.
        
        Args:
            dataset: DataFrame with protein sequences
            output_path: Path to save embeddings
            model_name: ESM model name
            layer: Layer to extract embeddings from
            pooling: Pooling method
            
        Returns:
            DataFrame with ESM embeddings
        """
        self.logger.info(f"Generating ESM embeddings using {model_name}...")
        
        try:
            import esm
            
            # Load ESM model
            model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
            model = model.to(self.device)
            model.eval()
            
            embeddings_data = []
            batch_converter = alphabet.get_batch_converter()
            
            # Process in batches
            batch_size = 32
            for i in range(0, len(dataset), batch_size):
                batch_df = dataset.iloc[i:i+batch_size]
                
                # Prepare batch
                batch_data = [(f"seq_{j}", seq) for j, seq in enumerate(batch_df['seq'])]
                batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
                batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
                
                batch_tokens = batch_tokens.to(self.device)
                
                # Generate embeddings
                with torch.no_grad():
                    results = model(batch_tokens, repr_layers=[layer], return_contacts=False)
                    token_representations = results["representations"][layer]
                
                # Apply pooling
                for j, tokens_len in enumerate(batch_lens):
                    if pooling == 'mean':
                        embedding = token_representations[j, 1:tokens_len-1].mean(0)
                    elif pooling == 'cls':
                        embedding = token_representations[j, 0]
                    else:
                        embedding = token_representations[j, tokens_len-2]
                    
                    embeddings_data.append({
                        'uniprot_id': batch_df.iloc[j]['uniprot_id'],
                        'esm_embedding': embedding.cpu().numpy()
                    })
                
                self.logger.info(f"Processed batch {i//batch_size + 1}/{(len(dataset) + batch_size - 1)//batch_size}")
            
            embeddings_df = pd.DataFrame(embeddings_data)
            save_dataframe(embeddings_df, output_path, 'feather')
            
            self.logger.info(f"ESM embeddings generated: {len(embeddings_df)} proteins")
            return embeddings_df
            
        except ImportError:
            self.logger.error("ESM not installed. Please install with: pip install fair-esm")
            raise
    
    def generate_t5_embeddings(self,
                              dataset: pd.DataFrame,
                              output_path: str,
                              model_name: str = "Rostlab/prot_t5_xl_half_uniref50-enc") -> pd.DataFrame:
        """
        Generate T5 embeddings for protein sequences.
        
        Args:
            dataset: DataFrame with protein sequences
            output_path: Path to save embeddings
            model_name: T5 model name
            
        Returns:
            DataFrame with T5 embeddings
        """
        self.logger.info(f"Generating T5 embeddings using {model_name}...")
        
        try:
            from transformers import T5Tokenizer, T5EncoderModel
            
            # Load T5 model
            tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
            model = T5EncoderModel.from_pretrained(model_name)
            model = model.to(self.device)
            model.eval()
            
            embeddings_data = []
            
            # Process sequences
            for _, row in dataset.iterrows():
                sequence = row['seq']
                uniprot_id = row['uniprot_id']
                
                # Tokenize sequence
                sequence_Example = ' '.join(list(sequence))
                encoded_input = tokenizer(sequence_Example, return_tensors='pt')
                input_ids = encoded_input['input_ids'].to(self.device)
                attention_mask = encoded_input['attention_mask'].to(self.device)
                
                # Generate embeddings
                with torch.no_grad():
                    output = model(input_ids=input_ids, attention_mask=attention_mask)
                    embedding = output.last_hidden_state.mean(dim=1).squeeze()
                
                embeddings_data.append({
                    'uniprot_id': uniprot_id,
                    't5_embedding': embedding.cpu().numpy()
                })
            
            embeddings_df = pd.DataFrame(embeddings_data)
            save_dataframe(embeddings_df, output_path, 'feather')
            
            self.logger.info(f"T5 embeddings generated: {len(embeddings_df)} proteins")
            return embeddings_df
            
        except ImportError:
            self.logger.error("Transformers not installed. Please install with: pip install transformers")
            raise
    
    def generate_unirep_embeddings(self,
                                  dataset: pd.DataFrame,
                                  output_path: str) -> pd.DataFrame:
        """
        Generate UniRep embeddings for protein sequences.
        
        Args:
            dataset: DataFrame with protein sequences
            output_path: Path to save embeddings
            
        Returns:
            DataFrame with UniRep embeddings
        """
        self.logger.info("Generating UniRep embeddings...")
        
        try:
            import unirep
            
            # Load UniRep model
            model = unirep.load()
            
            embeddings_data = []
            
            # Process sequences
            for _, row in dataset.iterrows():
                sequence = row['seq']
                uniprot_id = row['uniprot_id']
                
                # Generate embedding
                embedding = model.get_rep(sequence)
                
                embeddings_data.append({
                    'uniprot_id': uniprot_id,
                    'unirep_embedding': embedding
                })
            
            embeddings_df = pd.DataFrame(embeddings_data)
            save_dataframe(embeddings_df, output_path, 'feather')
            
            self.logger.info(f"UniRep embeddings generated: {len(embeddings_df)} proteins")
            return embeddings_df
            
        except ImportError:
            self.logger.error("UniRep not installed. Please install with: pip install unirep")
            raise
    
    def combine_embeddings(self,
                          embedding_files: Dict[str, str],
                          output_path: str) -> pd.DataFrame:
        """
        Combine multiple embedding types.
        
        Args:
            embedding_files: Dictionary mapping embedding type to file path
            output_path: Path to save combined embeddings
            
        Returns:
            DataFrame with combined embeddings
        """
        self.logger.info("Combining embeddings...")
        
        # Load all embedding files
        embeddings = {}
        for embedding_type, file_path in embedding_files.items():
            embeddings[embedding_type] = load_dataframe(file_path)
        
        # Find common proteins
        common_proteins = set.intersection(*[
            set(df['uniprot_id']) for df in embeddings.values()
        ])
        
        # Combine embeddings
        combined_data = []
        for uniprot_id in common_proteins:
            combined_embedding = {}
            combined_embedding['uniprot_id'] = uniprot_id
            
            for embedding_type, df in embeddings.items():
                protein_data = df[df['uniprot_id'] == uniprot_id].iloc[0]
                embedding_key = f"{embedding_type}_embedding"
                combined_embedding[embedding_key] = protein_data[embedding_key]
            
            combined_data.append(combined_embedding)
        
        combined_df = pd.DataFrame(combined_data)
        save_dataframe(combined_df, output_path, 'feather')
        
        self.logger.info(f"Embeddings combined: {len(combined_df)} proteins")
        return combined_df
    
    def normalize_embeddings(self,
                           embeddings_df: pd.DataFrame,
                           embedding_columns: List[str]) -> pd.DataFrame:
        """
        Normalize embeddings.
        
        Args:
            embeddings_df: DataFrame with embeddings
            embedding_columns: List of embedding column names
            
        Returns:
            DataFrame with normalized embeddings
        """
        self.logger.info("Normalizing embeddings...")
        
        normalized_df = embeddings_df.copy()
        
        for col in embedding_columns:
            if col in normalized_df.columns:
                embeddings = np.stack(normalized_df[col].values)
                
                # L2 normalization
                norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
                normalized_embeddings = embeddings / norms
                
                normalized_df[col] = list(normalized_embeddings)
        
        self.logger.info("Embeddings normalized")
        return normalized_df
    
    def reduce_embeddings_dimension(self,
                                  embeddings_df: pd.DataFrame,
                                  embedding_columns: List[str],
                                  target_dim: int = 512,
                                  method: str = 'pca') -> pd.DataFrame:
        """
        Reduce embedding dimensions.
        
        Args:
            embeddings_df: DataFrame with embeddings
            embedding_columns: List of embedding column names
            target_dim: Target dimension
            method: Dimensionality reduction method
            
        Returns:
            DataFrame with reduced embeddings
        """
        self.logger.info(f"Reducing embedding dimensions using {method}...")
        
        try:
            from sklearn.decomposition import PCA
            from sklearn.manifold import TSNE
            
            reduced_df = embeddings_df.copy()
            
            for col in embedding_columns:
                if col in reduced_df.columns:
                    embeddings = np.stack(reduced_df[col].values)
                    
                    if method == 'pca':
                        reducer = PCA(n_components=target_dim)
                    elif method == 'tsne':
                        reducer = TSNE(n_components=target_dim, random_state=42)
                    else:
                        raise ValueError(f"Unknown reduction method: {method}")
                    
                    reduced_embeddings = reducer.fit_transform(embeddings)
                    reduced_df[col] = list(reduced_embeddings)
            
            self.logger.info(f"Embeddings reduced to {target_dim} dimensions")
            return reduced_df
            
        except ImportError:
            self.logger.error("Scikit-learn not installed. Please install with: pip install scikit-learn")
            raise

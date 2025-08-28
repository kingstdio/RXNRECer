"""
ESM embedding model for RXNRECer
"""

import torch
import torch.nn as nn
import esm
from typing import List
import logging


class ESMEmbedding(nn.Module):
    """ESM embedding model with layer freezing."""
    
    def __init__(self, device: torch.device, freeze_layers: int = 32):
        """
        Initialize ESM embedding model.
        
        Args:
            device: Device to run model on
            freeze_layers: Number of layers to freeze
        """
        super(ESMEmbedding, self).__init__()
        
        self.logger = logging.getLogger(__name__)
        self.device = device
        self.freeze_layers = freeze_layers
        
        # Load ESM model
        self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        
        # Freeze specified layers
        self._freeze_layers()
        
        # Move to device
        self.esm_model = self.esm_model.to(self.device)
        self.esm_model.eval()
        
        self.logger.info(f"ESM model loaded with {freeze_layers} frozen layers")
    
    def _freeze_layers(self):
        """Freeze specified number of layers."""
        for name, param in self.esm_model.named_parameters():
            # Freeze layers up to freeze_layers
            if any(f"layer.{i}." in name for i in range(self.freeze_layers)):
                param.requires_grad = False
    
    def forward(self, sequences: List[str]) -> torch.Tensor:
        """
        Forward pass through ESM model.
        
        Args:
            sequences: List of protein sequences
            
        Returns:
            Tensor of sequence representations
        """
        # Prepare batch
        batch_converter = self.esm_alphabet.get_batch_converter()
        
        # Create sequence IDs
        seq_ids = [f"seq_{i}" for i in range(len(sequences))]
        batch_data = list(zip(seq_ids, sequences))
        
        # Convert to batch
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)
        
        # Move to device
        batch_tokens = batch_tokens.to(self.device)
        
        # Get embeddings
        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, 
                repr_layers=[33], 
                return_contacts=False
            )
        
        # Extract representations
        token_representations = results["representations"][33]
        
        # Pool to sequence representations
        sequence_representations = []
        for i, tokens_len in enumerate(batch_lens):
            # Take mean over sequence length (excluding special tokens)
            seq_repr = token_representations[i, 1:tokens_len-1].mean(0)
            sequence_representations.append(seq_repr)
        
        # Stack sequences
        sequence_representations = torch.stack(sequence_representations, dim=0)
        
        return sequence_representations.to(self.device)
    
    def get_embeddings(self, sequences: List[str], 
                      layer: int = 33,
                      pooling: str = 'mean') -> torch.Tensor:
        """
        Get embeddings from specific layer with specified pooling.
        
        Args:
            sequences: List of protein sequences
            layer: Layer to extract embeddings from
            pooling: Pooling method ('mean', 'cls', 'last')
            
        Returns:
            Tensor of embeddings
        """
        # Prepare batch
        batch_converter = self.esm_alphabet.get_batch_converter()
        seq_ids = [f"seq_{i}" for i in range(len(sequences))]
        batch_data = list(zip(seq_ids, sequences))
        
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)
        
        batch_tokens = batch_tokens.to(self.device)
        
        # Get embeddings
        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, 
                repr_layers=[layer], 
                return_contacts=False
            )
        
        token_representations = results["representations"][layer]
        
        # Apply pooling
        if pooling == 'mean':
            sequence_representations = []
            for i, tokens_len in enumerate(batch_lens):
                seq_repr = token_representations[i, 1:tokens_len-1].mean(0)
                sequence_representations.append(seq_repr)
            sequence_representations = torch.stack(sequence_representations, dim=0)
        
        elif pooling == 'cls':
            # Use CLS token representation
            sequence_representations = token_representations[:, 0]
        
        elif pooling == 'last':
            # Use last token representation
            sequence_representations = []
            for i, tokens_len in enumerate(batch_lens):
                seq_repr = token_representations[i, tokens_len-2]  # Last non-padding token
                sequence_representations.append(seq_repr)
            sequence_representations = torch.stack(sequence_representations, dim=0)
        
        else:
            raise ValueError(f"Unknown pooling method: {pooling}")
        
        return sequence_representations.to(self.device)
    
    def get_attention_weights(self, sequences: List[str], 
                            layer: int = 33) -> torch.Tensor:
        """
        Get attention weights from specified layer.
        
        Args:
            sequences: List of protein sequences
            layer: Layer to extract attention weights from
            
        Returns:
            Tensor of attention weights
        """
        # Prepare batch
        batch_converter = self.esm_alphabet.get_batch_converter()
        seq_ids = [f"seq_{i}" for i in range(len(sequences))]
        batch_data = list(zip(seq_ids, sequences))
        
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(self.device)
        
        # Get attention weights
        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, 
                repr_layers=[layer], 
                return_contacts=True
            )
        
        # Extract attention weights
        attention_weights = results["attentions"][layer]
        
        return attention_weights.to(self.device)
    
    def get_contact_map(self, sequences: List[str]) -> torch.Tensor:
        """
        Get contact map predictions.
        
        Args:
            sequences: List of protein sequences
            
        Returns:
            Tensor of contact maps
        """
        # Prepare batch
        batch_converter = self.esm_alphabet.get_batch_converter()
        seq_ids = [f"seq_{i}" for i in range(len(sequences))]
        batch_data = list(zip(seq_ids, sequences))
        
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(self.device)
        
        # Get contact maps
        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, 
                repr_layers=[33], 
                return_contacts=True
            )
        
        # Extract contact maps
        contact_maps = results["contacts"]
        
        return contact_maps.to(self.device)
    
    def unfreeze_layers(self, num_layers: int):
        """
        Unfreeze specified number of layers.
        
        Args:
            num_layers: Number of layers to unfreeze
        """
        for name, param in self.esm_model.named_parameters():
            if any(f"layer.{i}." in name for i in range(num_layers)):
                param.requires_grad = True
        
        self.logger.info(f"Unfroze {num_layers} layers")
    
    def freeze_all_layers(self):
        """Freeze all layers."""
        for param in self.esm_model.parameters():
            param.requires_grad = False
        
        self.logger.info("Froze all layers")
    
    def get_embedding_dim(self) -> int:
        """Get embedding dimension."""
        return self.esm_model.args.embed_dim

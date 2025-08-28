"""
BGRU model for RXNRECer
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional


class SelfAttention(nn.Module):
    """Self-attention mechanism."""
    
    def __init__(self, attention_size: int, input_dimensions: int):
        """
        Initialize self-attention layer.
        
        Args:
            attention_size: Size of attention mechanism
            input_dimensions: Input dimensions
        """
        super(SelfAttention, self).__init__()
        self.attention_size = attention_size
        
        # Attention parameters
        self.W = nn.Parameter(torch.Tensor(input_dimensions, self.attention_size))
        self.b = nn.Parameter(torch.Tensor(1, self.attention_size))
        self.u = nn.Parameter(torch.Tensor(self.attention_size, 1))
        
        self.reset_parameters()
    
    def reset_parameters(self):
        """Initialize parameters."""
        nn.init.xavier_normal_(self.W)
        nn.init.zeros_(self.b)
        nn.init.xavier_normal_(self.u)
    
    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through attention layer.
        
        Args:
            x: Input tensor [batch_size, seq_len, input_dim]
            mask: Attention mask [batch_size, seq_len]
            
        Returns:
            Attended output [batch_size, input_dim]
        """
        # Calculate attention scores
        # et: [batch_size, seq_len, attention_size]
        et = torch.tanh(torch.matmul(x, self.W) + self.b)
        
        # at: [batch_size, seq_len, 1]
        at = torch.matmul(et, self.u)
        
        # Apply mask if provided
        if mask is not None:
            at = at.masked_fill(mask.unsqueeze(-1) == 0, float('-inf'))
        
        # Softmax to get attention weights
        # at: [batch_size, seq_len, 1]
        at = F.softmax(at, dim=1)
        
        # Apply attention weights
        # output: [batch_size, input_dim]
        output = torch.sum(x * at, dim=1)
        
        return output


class BGRUModel(nn.Module):
    """Bidirectional GRU model with attention."""
    
    def __init__(self, input_dimensions: int, gru_h_size: int = 512,
                 attention_size: int = 32, dropout: float = 0.2,
                 output_dimensions: int = 10479, device: torch.device = None):
        """
        Initialize BGRU model.
        
        Args:
            input_dimensions: Input feature dimensions
            gru_h_size: GRU hidden size
            attention_size: Attention mechanism size
            dropout: Dropout rate
            output_dimensions: Output dimensions (number of reactions)
            device: Device to run model on
        """
        super(BGRUModel, self).__init__()
        
        self.input_dimensions = input_dimensions
        self.gru_h_size = gru_h_size
        self.attention_size = attention_size
        self.dropout = dropout
        self.output_dimensions = output_dimensions
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        # Bidirectional GRU
        self.gru = nn.GRU(
            input_size=input_dimensions,
            hidden_size=gru_h_size,
            num_layers=1,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if 1 > 1 else 0
        )
        
        # Attention mechanism
        self.attention = SelfAttention(
            attention_size=attention_size,
            input_dimensions=gru_h_size * 2  # Bidirectional
        )
        
        # Output layer
        self.output_layer = nn.Linear(gru_h_size * 2, output_dimensions)
        
        # Dropout layer
        self.dropout_layer = nn.Dropout(dropout)
        
        # Move to device
        self.to(self.device)
    
    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Forward pass through BGRU model.
        
        Args:
            x: Input tensor [batch_size, seq_len, input_dim]
            mask: Sequence mask [batch_size, seq_len]
            
        Returns:
            Output logits [batch_size, output_dimensions]
        """
        # GRU forward pass
        # gru_out: [batch_size, seq_len, gru_h_size * 2]
        gru_out, _ = self.gru(x)
        
        # Apply attention
        # attended: [batch_size, gru_h_size * 2]
        attended = self.attention(gru_out, mask)
        
        # Apply dropout
        attended = self.dropout_layer(attended)
        
        # Output layer
        # output: [batch_size, output_dimensions]
        output = self.output_layer(attended)
        
        return output
    
    def get_attention_weights(self, x: torch.Tensor, 
                            mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Get attention weights for visualization.
        
        Args:
            x: Input tensor
            mask: Sequence mask
            
        Returns:
            Attention weights [batch_size, seq_len]
        """
        # GRU forward pass
        gru_out, _ = self.gru(x)
        
        # Calculate attention scores
        et = torch.tanh(torch.matmul(gru_out, self.attention.W) + self.attention.b)
        at = torch.matmul(et, self.attention.u).squeeze(-1)
        
        # Apply mask if provided
        if mask is not None:
            at = at.masked_fill(mask == 0, float('-inf'))
        
        # Softmax to get attention weights
        attention_weights = F.softmax(at, dim=1)
        
        return attention_weights
    
    def get_embeddings(self, x: torch.Tensor, 
                      mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Get intermediate embeddings.
        
        Args:
            x: Input tensor
            mask: Sequence mask
            
        Returns:
            Embeddings [batch_size, gru_h_size * 2]
        """
        # GRU forward pass
        gru_out, _ = self.gru(x)
        
        # Apply attention
        attended = self.attention(gru_out, mask)
        
        return attended


class BGRUWithESM(nn.Module):
    """Combined ESM + BGRU model."""
    
    def __init__(self, esm_embedding: nn.Module, bgru_model: BGRUModel):
        """
        Initialize combined model.
        
        Args:
            esm_embedding: ESM embedding model
            bgru_model: BGRU model
        """
        super(BGRUWithESM, self).__init__()
        
        self.esm_embedding = esm_embedding
        self.bgru_model = bgru_model
    
    def forward(self, sequences: list) -> torch.Tensor:
        """
        Forward pass through combined model.
        
        Args:
            sequences: List of protein sequences
            
        Returns:
            Output logits
        """
        # Get ESM embeddings
        esm_embeddings = self.esm_embedding(sequences)
        
        # Reshape for GRU (add sequence dimension)
        # esm_embeddings: [batch_size, esm_dim] -> [batch_size, 1, esm_dim]
        esm_embeddings = esm_embeddings.unsqueeze(1)
        
        # Pass through BGRU
        output = self.bgru_model(esm_embeddings)
        
        return output
    
    def get_attention_weights(self, sequences: list) -> torch.Tensor:
        """Get attention weights."""
        esm_embeddings = self.esm_embedding(sequences)
        esm_embeddings = esm_embeddings.unsqueeze(1)
        return self.bgru_model.get_attention_weights(esm_embeddings)
    
    def get_embeddings(self, sequences: list) -> torch.Tensor:
        """Get intermediate embeddings."""
        esm_embeddings = self.esm_embedding(sequences)
        esm_embeddings = esm_embeddings.unsqueeze(1)
        return self.bgru_model.get_embeddings(esm_embeddings)


def create_bgru_model(input_dimensions: int, 
                     gru_h_size: int = 512,
                     attention_size: int = 32,
                     dropout: float = 0.2,
                     output_dimensions: int = 10479,
                     device: Optional[torch.device] = None) -> BGRUModel:
    """
    Create BGRU model with default parameters.
    
    Args:
        input_dimensions: Input feature dimensions
        gru_h_size: GRU hidden size
        attention_size: Attention mechanism size
        dropout: Dropout rate
        output_dimensions: Output dimensions
        device: Device to run model on
        
    Returns:
        Initialized BGRU model
    """
    return BGRUModel(
        input_dimensions=input_dimensions,
        gru_h_size=gru_h_size,
        attention_size=attention_size,
        dropout=dropout,
        output_dimensions=output_dimensions,
        device=device
    )

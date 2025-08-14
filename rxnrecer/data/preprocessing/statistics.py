"""
Dataset statistics and analysis utilities for RXNRECer
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import logging

from ...utils.bio_utils import get_protein_properties


class DatasetStatistics:
    """Dataset statistics and analysis utilities."""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
    
    def analyze_sequence_statistics(self, dataset: pd.DataFrame) -> Dict:
        """
        Analyze sequence statistics.
        
        Args:
            dataset: Dataset with protein sequences
            
        Returns:
            Dictionary with sequence statistics
        """
        self.logger.info("Analyzing sequence statistics...")
        
        sequences = dataset['seq'].tolist()
        
        # Basic statistics
        lengths = [len(seq) for seq in sequences]
        stats = {
            'total_sequences': len(sequences),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'std_length': np.std(lengths),
            'length_distribution': {
                'short': len([l for l in lengths if l < 100]),
                'medium': len([l for l in lengths if 100 <= l < 500]),
                'long': len([l for l in lengths if l >= 500])
            }
        }
        
        # Amino acid composition
        aa_composition = self._calculate_aa_composition(sequences)
        stats['amino_acid_composition'] = aa_composition
        
        # Protein properties
        properties = self._calculate_protein_properties(sequences)
        stats['protein_properties'] = properties
        
        return stats
    
    def _calculate_aa_composition(self, sequences: List[str]) -> Dict[str, float]:
        """Calculate amino acid composition."""
        aa_counts = {}
        total_aa = 0
        
        for seq in sequences:
            for aa in seq:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
                total_aa += 1
        
        # Convert to percentages
        aa_composition = {aa: count/total_aa*100 for aa, count in aa_counts.items()}
        return aa_composition
    
    def _calculate_protein_properties(self, sequences: List[str]) -> Dict:
        """Calculate protein properties statistics."""
        properties_list = [get_protein_properties(seq) for seq in sequences]
        
        # Extract properties
        mw_list = [p['molecular_weight'] for p in properties_list]
        pi_list = [p['isoelectric_point'] for p in properties_list]
        hydrophobicity_list = [p['hydrophobicity'] for p in properties_list]
        charge_list = [p['charge'] for p in properties_list]
        
        return {
            'molecular_weight': {
                'mean': np.mean(mw_list),
                'std': np.std(mw_list),
                'min': min(mw_list),
                'max': max(mw_list)
            },
            'isoelectric_point': {
                'mean': np.mean(pi_list),
                'std': np.std(pi_list),
                'min': min(pi_list),
                'max': max(pi_list)
            },
            'hydrophobicity': {
                'mean': np.mean(hydrophobicity_list),
                'std': np.std(hydrophobicity_list),
                'min': min(hydrophobicity_list),
                'max': max(hydrophobicity_list)
            },
            'charge': {
                'mean': np.mean(charge_list),
                'std': np.std(charge_list),
                'min': min(charge_list),
                'max': max(charge_list)
            }
        }
    
    def analyze_reaction_statistics(self, dataset: pd.DataFrame) -> Dict:
        """
        Analyze reaction statistics.
        
        Args:
            dataset: Dataset with reaction information
            
        Returns:
            Dictionary with reaction statistics
        """
        self.logger.info("Analyzing reaction statistics...")
        
        # Reaction frequency
        reaction_counts = dataset['reaction_id'].value_counts()
        
        stats = {
            'total_reactions': len(dataset),
            'unique_reactions': len(reaction_counts),
            'reaction_frequency': {
                'min': reaction_counts.min(),
                'max': reaction_counts.max(),
                'mean': reaction_counts.mean(),
                'median': reaction_counts.median()
            },
            'top_reactions': reaction_counts.head(10).to_dict(),
            'reaction_distribution': {
                'single_occurrence': (reaction_counts == 1).sum(),
                'rare': ((reaction_counts > 1) & (reaction_counts <= 5)).sum(),
                'common': ((reaction_counts > 5) & (reaction_counts <= 50)).sum(),
                'frequent': (reaction_counts > 50).sum()
            }
        }
        
        return stats
    
    def analyze_dataset_balance(self, dataset: pd.DataFrame) -> Dict:
        """
        Analyze dataset balance.
        
        Args:
            dataset: Dataset to analyze
            
        Returns:
            Dictionary with balance statistics
        """
        self.logger.info("Analyzing dataset balance...")
        
        # Check for class imbalance
        if 'reaction_id' in dataset.columns:
            reaction_counts = dataset['reaction_id'].value_counts()
            
            balance_stats = {
                'total_samples': len(dataset),
                'total_classes': len(reaction_counts),
                'class_imbalance_ratio': reaction_counts.max() / reaction_counts.min(),
                'gini_coefficient': self._calculate_gini_coefficient(reaction_counts.values),
                'entropy': self._calculate_entropy(reaction_counts.values)
            }
        else:
            balance_stats = {'error': 'No reaction_id column found'}
        
        return balance_stats
    
    def _calculate_gini_coefficient(self, values: np.ndarray) -> float:
        """Calculate Gini coefficient for class imbalance."""
        sorted_values = np.sort(values)
        n = len(sorted_values)
        cumsum = np.cumsum(sorted_values)
        return (n + 1 - 2 * np.sum(cumsum) / cumsum[-1]) / n
    
    def _calculate_entropy(self, values: np.ndarray) -> float:
        """Calculate entropy for class distribution."""
        probabilities = values / np.sum(values)
        return -np.sum(probabilities * np.log2(probabilities))
    
    def generate_statistics_report(self, 
                                 dataset: pd.DataFrame,
                                 output_path: str) -> Dict:
        """
        Generate comprehensive statistics report.
        
        Args:
            dataset: Dataset to analyze
            output_path: Path to save report
            
        Returns:
            Complete statistics report
        """
        self.logger.info("Generating statistics report...")
        
        report = {
            'dataset_info': {
                'total_samples': len(dataset),
                'columns': list(dataset.columns),
                'missing_values': dataset.isnull().sum().to_dict()
            },
            'sequence_statistics': self.analyze_sequence_statistics(dataset),
            'reaction_statistics': self.analyze_reaction_statistics(dataset),
            'balance_statistics': self.analyze_dataset_balance(dataset)
        }
        
        # Save report
        import json
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        self.logger.info(f"Statistics report saved to {output_path}")
        return report
    
    def plot_sequence_length_distribution(self, 
                                        dataset: pd.DataFrame,
                                        output_path: str):
        """Plot sequence length distribution."""
        lengths = dataset['seq'].str.len()
        
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Sequence Length')
        plt.ylabel('Frequency')
        plt.title('Protein Sequence Length Distribution')
        plt.grid(True, alpha=0.3)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Length distribution plot saved to {output_path}")
    
    def plot_reaction_frequency(self, 
                              dataset: pd.DataFrame,
                              output_path: str,
                              top_n: int = 20):
        """Plot reaction frequency distribution."""
        reaction_counts = dataset['reaction_id'].value_counts().head(top_n)
        
        plt.figure(figsize=(12, 8))
        reaction_counts.plot(kind='bar')
        plt.xlabel('Reaction ID')
        plt.ylabel('Frequency')
        plt.title(f'Top {top_n} Most Frequent Reactions')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Reaction frequency plot saved to {output_path}")

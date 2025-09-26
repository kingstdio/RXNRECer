"""
Evaluation utility functions for RXNRECer
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Union
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, average_precision_score, confusion_matrix,
    classification_report
)
from sklearn.model_selection import train_test_split
import math


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray, 
                     y_prob: Optional[np.ndarray] = None,
                     average: str = 'weighted') -> Dict[str, float]:
    """
    Calculate comprehensive evaluation metrics.
    
    Args:
        y_true: True labels
        y_pred: Predicted labels
        y_prob: Predicted probabilities (optional)
        average: Averaging method for multi-class metrics
        
    Returns:
        Dictionary with evaluation metrics
    """
    metrics = {}
    
    # Basic classification metrics
    metrics['accuracy'] = accuracy_score(y_true, y_pred)
    metrics['precision'] = precision_score(y_true, y_pred, average=average, zero_division=0)
    metrics['recall'] = recall_score(y_true, y_pred, average=average, zero_division=0)
    metrics['f1'] = f1_score(y_true, y_pred, average=average, zero_division=0)
    
    # Probability-based metrics
    if y_prob is not None:
        try:
            if len(y_prob.shape) == 1:
                # Binary classification
                metrics['roc_auc'] = roc_auc_score(y_true, y_prob)
                metrics['pr_auc'] = average_precision_score(y_true, y_prob)
            else:
                # Multi-class classification
                metrics['roc_auc'] = roc_auc_score(y_true, y_prob, multi_class='ovr', average=average)
                metrics['pr_auc'] = average_precision_score(y_true, y_prob, average=average)
        except Exception as e:
            print(f"Warning: Could not calculate probability-based metrics: {e}")
    
    return metrics


def calculate_reaction_metrics(ground_truth: pd.DataFrame, predictions: pd.DataFrame,
                             reaction_col: str = 'reaction_id',
                             top_k: int = 5) -> Dict[str, float]:
    """
    Calculate reaction prediction specific metrics.
    
    Args:
        ground_truth: Ground truth DataFrame
        predictions: Predictions DataFrame
        reaction_col: Column name for reaction IDs
        top_k: Top-k predictions to consider
        
    Returns:
        Dictionary with reaction-specific metrics
    """
    metrics = {}
    
    # Ensure both DataFrames have the same index
    common_idx = ground_truth.index.intersection(predictions.index)
    gt = ground_truth.loc[common_idx]
    pred = predictions.loc[common_idx]
    
    # Top-k accuracy
    if f'top_{top_k}_predictions' in pred.columns:
        top_k_correct = 0
        for idx in common_idx:
            true_rxn = gt.loc[idx, reaction_col]
            pred_rxns = pred.loc[idx, f'top_{top_k}_predictions']
            if isinstance(pred_rxns, str):
                pred_rxns = pred_rxns.split(';')
            if true_rxn in pred_rxns:
                top_k_correct += 1
        
        metrics[f'top_{top_k}_accuracy'] = top_k_correct / len(common_idx)
    
    # Mean Reciprocal Rank (MRR)
    if 'prediction_rank' in pred.columns:
        reciprocal_ranks = []
        for idx in common_idx:
            rank = pred.loc[idx, 'prediction_rank']
            if rank > 0:
                reciprocal_ranks.append(1.0 / rank)
            else:
                reciprocal_ranks.append(0.0)
        
        metrics['mrr'] = np.mean(reciprocal_ranks)
    
    return metrics


def split_dataset(dataset: pd.DataFrame, test_ratio: float = 0.2,
                 split_method: str = 'random', 
                 stratify_col: Optional[str] = None,
                 random_state: int = 42) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split dataset into train and test sets.
    
    Args:
        dataset: Input dataset
        test_ratio: Ratio of test set
        split_method: Split method ('random', 'protein', 'reaction')
        stratify_col: Column to stratify on
        random_state: Random seed
        
    Returns:
        Tuple of (train_df, test_df)
    """
    if split_method == 'random':
        if stratify_col:
            train_df, test_df = train_test_split(
                dataset, test_size=test_ratio, 
                stratify=dataset[stratify_col],
                random_state=random_state
            )
        else:
            train_df, test_df = train_test_split(
                dataset, test_size=test_ratio,
                random_state=random_state
            )
    
    elif split_method == 'protein':
        # Split by unique proteins
        proteins = dataset['uniprot_id'].unique()
        test_proteins = np.random.choice(
            proteins, 
            size=int(len(proteins) * test_ratio),
            replace=False
        )
        
        test_df = dataset[dataset['uniprot_id'].isin(test_proteins)]
        train_df = dataset[~dataset['uniprot_id'].isin(test_proteins)]
    
    elif split_method == 'reaction':
        # Split by unique reactions
        reactions = dataset['reaction_id'].unique()
        test_reactions = np.random.choice(
            reactions,
            size=int(len(reactions) * test_ratio),
            replace=False
        )
        
        test_df = dataset[dataset['reaction_id'].isin(test_reactions)]
        train_df = dataset[~dataset['reaction_id'].isin(test_reactions)]
    
    else:
        raise ValueError(f"Unknown split method: {split_method}")
    
    # Add split labels
    train_df = train_df.copy()
    test_df = test_df.copy()
    train_df['split'] = 'train'
    test_df['split'] = 'test'
    
    return train_df, test_df


def calculate_confidence_intervals(metrics: Dict[str, List[float]], 
                                 confidence_level: float = 0.95) -> Dict[str, Dict[str, float]]:
    """
    Calculate confidence intervals for metrics.
    
    Args:
        metrics: Dictionary with metric names and lists of values
        confidence_level: Confidence level (0.95 for 95% CI)
        
    Returns:
        Dictionary with confidence intervals
    """
    intervals = {}
    
    for metric_name, values in metrics.items():
        if len(values) < 2:
            continue
            
        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)
        n = len(values)
        
        # t-distribution critical value
        t_critical = 1.96  # Approximate for 95% CI with large n
        
        margin_of_error = t_critical * (std_val / np.sqrt(n))
        
        intervals[metric_name] = {
            'mean': mean_val,
            'std': std_val,
            'lower': mean_val - margin_of_error,
            'upper': mean_val + margin_of_error,
            'n': n
        }
    
    return intervals


def format_classification_report(y_true: np.ndarray, y_pred: np.ndarray,
                               target_names: Optional[List[str]] = None) -> str:
    """
    Format classification report.
    
    Args:
        y_true: True labels
        y_pred: Predicted labels
        target_names: Names for target classes
        
    Returns:
        Formatted classification report
    """
    report = classification_report(y_true, y_pred, target_names=target_names)
    return report


def calculate_hit_rate(predictions: pd.DataFrame, ground_truth: pd.DataFrame,
                      top_k: int = 1, reaction_col: str = 'reaction_id') -> float:
    """
    Calculate hit rate (percentage of correct predictions in top-k).
    
    Args:
        predictions: Predictions DataFrame
        ground_truth: Ground truth DataFrame
        top_k: Top-k predictions to consider
        reaction_col: Column name for reaction IDs
        
    Returns:
        Hit rate
    """
    hits = 0
    total = 0
    
    for idx in ground_truth.index:
        if idx not in predictions.index:
            continue
            
        true_rxn = ground_truth.loc[idx, reaction_col]
        
        # Get top-k predictions
        if f'top_{top_k}_predictions' in predictions.columns:
            pred_rxns = predictions.loc[idx, f'top_{top_k}_predictions']
            if isinstance(pred_rxns, str):
                pred_rxns = pred_rxns.split(';')
            
            if true_rxn in pred_rxns[:top_k]:
                hits += 1
            total += 1
    
    return hits / total if total > 0 else 0.0


def calculate_coverage(predictions: pd.DataFrame, 
                      unique_reactions: Optional[List[str]] = None) -> float:
    """
    Calculate prediction coverage (percentage of unique reactions predicted).
    
    Args:
        predictions: Predictions DataFrame
        unique_reactions: List of all unique reactions (optional)
        
    Returns:
        Coverage percentage
    """
    if unique_reactions is None:
        # Extract from predictions
        all_pred_rxns = []
        for col in predictions.columns:
            if col.startswith('top_') and col.endswith('_predictions'):
                pred_rxns = predictions[col].dropna()
                for pred_str in pred_rxns:
                    if isinstance(pred_str, str):
                        all_pred_rxns.extend(pred_str.split(';'))
        
        unique_reactions = list(set(all_pred_rxns))
    
    predicted_reactions = set()
    for col in predictions.columns:
        if col.startswith('top_') and col.endswith('_predictions'):
            pred_rxns = predictions[col].dropna()
            for pred_str in pred_rxns:
                if isinstance(pred_str, str):
                    predicted_reactions.update(pred_str.split(';'))
    
    return len(predicted_reactions) / len(unique_reactions) if unique_reactions else 0.0


def evaluate_ensemble_predictions(predictions_list: List[pd.DataFrame],
                                ground_truth: pd.DataFrame,
                                ensemble_method: str = 'voting',
                                weights: Optional[List[float]] = None) -> Dict[str, float]:
    """
    Evaluate ensemble predictions.
    
    Args:
        predictions_list: List of prediction DataFrames
        ground_truth: Ground truth DataFrame
        ensemble_method: Ensemble method ('voting', 'weighted_voting', 'stacking')
        weights: Weights for weighted voting
        
    Returns:
        Dictionary with ensemble evaluation metrics
    """
    if ensemble_method == 'voting':
        # Simple majority voting
        ensemble_pred = voting_ensemble(predictions_list)
    elif ensemble_method == 'weighted_voting':
        # Weighted voting
        if weights is None:
            weights = [1.0] * len(predictions_list)
        ensemble_pred = weighted_voting_ensemble(predictions_list, weights)
    else:
        raise ValueError(f"Unknown ensemble method: {ensemble_method}")
    
    # Calculate metrics
    metrics = calculate_reaction_metrics(ground_truth, ensemble_pred)
    
    return metrics



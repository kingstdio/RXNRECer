"""
Biological utility functions for RXNRECer
"""

import os
import re
import subprocess
import tempfile
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML
import requests
from requests.exceptions import RequestException


def split_compound_str(compounds_str: str, remove_reaction_coefficient: bool = False) -> List[str]:
    """
    Split compound string into list of compounds.
    
    Args:
        compounds_str: String containing compounds
        remove_reaction_coefficient: Whether to remove reaction coefficients
        
    Returns:
        List of compound strings
    """
    if not compounds_str or compounds_str == '-':
        return []
    
    # Split by common separators
    compounds = re.split(r'[+\s]+', compounds_str.strip())
    
    if remove_reaction_coefficient:
        # Remove numeric coefficients
        compounds = [re.sub(r'^\d+\.?\d*\s*', '', comp) for comp in compounds]
    
    # Filter out empty strings
    compounds = [comp.strip() for comp in compounds if comp.strip()]
    
    return compounds


def get_blast_results(train_df: pd.DataFrame, test_df: pd.DataFrame, 
                     k: int = 1, evalue: float = 1e-5) -> pd.DataFrame:
    """
    Perform BLAST search between train and test datasets.
    
    Args:
        train_df: Training dataset DataFrame
        test_df: Test dataset DataFrame
        k: Number of top hits to return
        evalue: E-value threshold
        
    Returns:
        DataFrame with BLAST results
    """
    with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_train, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as fasta_test, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.tsv') as res_blast, \
         tempfile.NamedTemporaryFile(delete=True, suffix='.dmnd') as db_dmnd:
        
        # Convert DataFrames to FASTA
        dataframe_to_fasta(train_df, fasta_train.name)
        dataframe_to_fasta(test_df, fasta_test.name)
        
        # Build DIAMOND database
        cmd1 = ["diamond", "makedb", "--in", fasta_train.name, "-d", db_dmnd.name, "--quiet"]
        subprocess.run(cmd1, check=True)
        
        # Run BLAST search
        cmd2 = ["diamond", "blastp", "-d", db_dmnd.name, "-q", fasta_test.name, 
                "-o", res_blast.name, "-b5", "-c1", "-k", str(k), 
                "-e", str(evalue), "--quiet"]
        subprocess.run(cmd2, check=True)
        
        # Read results
        columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        results = pd.read_csv(res_blast.name, sep='\t', names=columns)
        
    return results


def dataframe_to_fasta(df: pd.DataFrame, output_file: str, 
                      id_col: str = "uniprot_id", seq_col: str = "seq") -> None:
    """Convert DataFrame to FASTA file."""
    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row[id_col]}\n")
            f.write(f"{row[seq_col]}\n")


def calculate_sequence_similarity(seq1: str, seq2: str, method: str = 'identity') -> float:
    """
    Calculate sequence similarity between two sequences.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        method: Similarity method ('identity', 'blosum62', 'pam250')
        
    Returns:
        Similarity score
    """
    if method == 'identity':
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)
    else:
        # TODO: Implement other similarity methods
        raise NotImplementedError(f"Method {method} not implemented")


def get_protein_properties(sequence: str) -> Dict[str, float]:
    """
    Calculate basic protein properties.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Dictionary with protein properties
    """
    # Amino acid frequencies
    aa_freq = {}
    for aa in sequence:
        aa_freq[aa] = aa_freq.get(aa, 0) + 1
    
    # Normalize frequencies
    seq_len = len(sequence)
    aa_freq = {aa: count/seq_len for aa, count in aa_freq.items()}
    
    # Calculate properties
    properties = {
        'length': seq_len,
        'molecular_weight': calculate_molecular_weight(sequence),
        'isoelectric_point': calculate_isoelectric_point(sequence),
        'hydrophobicity': calculate_hydrophobicity(sequence),
        'charge': calculate_charge(sequence)
    }
    
    return properties


def calculate_molecular_weight(sequence: str) -> float:
    """Calculate molecular weight of protein sequence."""
    # Amino acid molecular weights (Da)
    aa_weights = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    
    weight = sum(aa_weights.get(aa, 0) for aa in sequence.upper())
    return weight - 18.0 * (len(sequence) - 1)  # Subtract water molecules


def calculate_isoelectric_point(sequence: str) -> float:
    """Calculate isoelectric point of protein sequence."""
    # pKa values for amino acids
    pka_values = {
        'D': 3.65, 'E': 4.25, 'H': 6.0, 'K': 10.53, 'R': 12.48,
        'Y': 10.07, 'C': 8.18
    }
    
    # Count charged residues
    charged_residues = {aa: sequence.upper().count(aa) for aa in pka_values}
    
    # Simple calculation (can be improved with more sophisticated methods)
    net_charge = (charged_residues['D'] + charged_residues['E'] - 
                  charged_residues['K'] - charged_residues['R'] - 
                  charged_residues['H'])
    
    # Approximate pI calculation
    if net_charge > 0:
        return 10.0
    elif net_charge < 0:
        return 4.0
    else:
        return 7.0


def calculate_hydrophobicity(sequence: str) -> float:
    """Calculate hydrophobicity of protein sequence."""
    # Kyte-Doolittle hydrophobicity scale
    hydrophobicity_scores = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    scores = [hydrophobicity_scores.get(aa, 0) for aa in sequence.upper()]
    return sum(scores) / len(scores)


def calculate_charge(sequence: str) -> float:
    """Calculate net charge of protein sequence at pH 7."""
    # Charge at pH 7
    charges = {
        'D': -1, 'E': -1, 'H': 0.1, 'K': 1, 'R': 1, 'Y': 0
    }
    
    net_charge = sum(charges.get(aa, 0) for aa in sequence.upper())
    return net_charge


def download_uniprot_entry(uniprot_id: str, save_path: Optional[str] = None) -> str:
    """
    Download UniProt entry.
    
    Args:
        uniprot_id: UniProt ID
        save_path: Path to save file (optional)
        
    Returns:
        Path to downloaded file
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    
    if save_path is None:
        save_path = f"{uniprot_id}.fasta"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        with open(save_path, 'w') as f:
            f.write(response.text)
        
        return save_path
    except RequestException as e:
        raise Exception(f"Failed to download UniProt entry {uniprot_id}: {e}")


def parse_fasta_header(header: str) -> Dict[str, str]:
    """
    Parse FASTA header to extract information.
    
    Args:
        header: FASTA header line
        
    Returns:
        Dictionary with parsed information
    """
    info = {}
    
    # Extract ID
    if '|' in header:
        parts = header.split('|')
        if len(parts) >= 2:
            info['database'] = parts[0].lstrip('>')
            info['id'] = parts[1]
            if len(parts) >= 3:
                info['entry_name'] = parts[2]
    
    # Extract description
    if ' ' in header:
        description = header.split(' ', 1)[1]
        info['description'] = description
    
    return info


def validate_sequence(sequence: str) -> bool:
    """
    Validate protein sequence.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        True if valid, False otherwise
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return all(aa in valid_aa for aa in sequence.upper())


def clean_sequence(sequence: str) -> str:
    """
    Clean protein sequence by removing invalid characters.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Cleaned sequence
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    return ''.join(aa for aa in sequence.upper() if aa in valid_aa)

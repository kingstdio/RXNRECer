"""
Biological utility functions for RXNRECer
"""

import sys
import os
import re
import subprocess
import tempfile
import pandas as pd
import numpy as np
import time
from typing import List, Dict, Tuple, Optional, Union
from collections import Counter
from pandarallel import pandarallel
import requests
from requests.exceptions import RequestException

from rxnrecer.config import config as cfg
from rxnrecer.utils import file_utils as ftool
from rxnrecer.lib.rxn.Reaction import Reaction

# Valid amino acids constant
VALID_AA = frozenset('ACDEFGHIKLMNPQRSTVWY')

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
    
    # Split by plus sign with optional surrounding whitespace
    compounds = re.split(r'\s*\+\s*', compounds_str.strip())
    
    if remove_reaction_coefficient:
        # Remove numeric coefficients
        compounds = [re.sub(r'^\d+\.?\d*\s*', '', comp) for comp in compounds]
    
    # Filter out empty strings
    return [comp.strip() for comp in compounds if comp.strip()]


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
        ftool.dataframe_to_fasta(train_df, fasta_train.name)
        ftool.dataframe_to_fasta(test_df, fasta_test.name)
        
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
    seq_len = len(sequence)
    if seq_len == 0:
        return {
            'length': 0, 'molecular_weight': 0.0, 'isoelectric_point': 7.0,
            'hydrophobicity': 0.0, 'charge': 0.0
        }
        
    aa_counts = Counter(sequence.upper())
    
    properties = {
        'length': seq_len,
        'molecular_weight': calculate_molecular_weight(sequence, aa_counts),
        'isoelectric_point': calculate_isoelectric_point(sequence, aa_counts),
        'hydrophobicity': calculate_hydrophobicity(sequence, aa_counts),
        'charge': calculate_charge(sequence, aa_counts)
    }
    
    return properties


def calculate_molecular_weight(sequence: str, aa_counts: Optional[Counter] = None) -> float:
    """Calculate molecular weight of protein sequence."""
    if not sequence:
        return 0.0
    aa_weights = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }
    counts = aa_counts if aa_counts is not None else Counter(sequence.upper())
    weight = sum(aa_weights.get(aa, 0) * count for aa, count in counts.items())
    return weight - 18.0 * (len(sequence) - 1)  # Subtract water molecules


def calculate_isoelectric_point(sequence: str, aa_counts: Optional[Counter] = None) -> float:
    """Calculate isoelectric point of protein sequence."""
    pka_values = {
        'D': 3.65, 'E': 4.25, 'H': 6.0, 'K': 10.53, 'R': 12.48,
        'Y': 10.07, 'C': 8.18
    }
    counts = aa_counts if aa_counts is not None else Counter(sequence.upper())
    
    # Simple calculation (can be improved with more sophisticated methods)
    net_charge = (counts['D'] + counts['E'] - counts['K'] - counts['R'] - counts['H'])
    
    # Approximate pI calculation
    if net_charge > 0:
        return 10.0
    elif net_charge < 0:
        return 4.0
    else:
        return 7.0


def calculate_hydrophobicity(sequence: str, aa_counts: Optional[Counter] = None) -> float:
    """Calculate hydrophobicity of protein sequence."""
    if not sequence:
        return 0.0
    hydrophobicity_scores = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    counts = aa_counts if aa_counts is not None else Counter(sequence.upper())
    total_score = sum(hydrophobicity_scores.get(aa, 0) * count for aa, count in counts.items())
    return total_score / len(sequence)


def calculate_charge(sequence: str, aa_counts: Optional[Counter] = None) -> float:
    """Calculate net charge of protein sequence at pH 7."""
    charges = {
        'D': -1, 'E': -1, 'H': 0.1, 'K': 1, 'R': 1, 'Y': 0
    }
    counts = aa_counts if aa_counts is not None else Counter(sequence.upper())
    net_charge = sum(charges.get(aa, 0) * count for aa, count in counts.items())
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
    return all(aa in VALID_AA for aa in sequence.upper())


def clean_sequence(sequence: str) -> str:
    """
    Clean protein sequence by removing invalid characters.
    
    Args:
        sequence: Protein sequence
        
    Returns:
        Cleaned sequence
    """
    return ''.join(aa for aa in sequence.upper() if aa in VALID_AA)


def get_rxn_detail(rxn_id: str, rxn_bank: pd.DataFrame) -> Optional[Reaction]:
    """
    Get detailed reaction information by reaction ID, returns Reaction object

    Args:
        rxn_id (str): Reaction ID string, e.g., 'RHEA:12345'
        rxn_bank (pd.DataFrame): Reaction database DataFrame, should contain columns:
            - reaction_id: Reaction ID
            - equation: Reaction equation
            - equation_chebi: ChEBI format reaction equation
            - equation_smiles: SMILES format reaction equation
            - ec_number: EC number

    Returns:
        Reaction: Reaction object containing complete reaction information
            Returns None if reaction ID is invalid or not found
    """
    if not rxn_id or rxn_id == '-':
        return None

    match = rxn_bank[rxn_bank.reaction_id == rxn_id.strip()]
    if match.empty:
        return None

    row = match.iloc[0]
    
    try:
        # Create Reaction object
        reaction = Reaction(
            rxn_smiles=row.equation_smiles,
            rxn_equation=row.equation,
            rxn_equation_ref_chebi=row.equation_chebi,
            rxn_id=row.reaction_id,
            rxn_ec=row.ec_number
        )
        return reaction
    except Exception as e:
        print(f"Warning: Unable to create reaction object {rxn_id}: {e}")
        return None


def get_rxn_details_from_rxn_json(rxn_ids: str) -> pd.DataFrame:
    """
    Get reaction details from JSON files with support for multiple separators and exception handling

    Args:
        rxn_ids (str): Reaction ID string supporting multiple separators (; | ,)

    Returns:
        pd.DataFrame: DataFrame containing reaction details
    """
    if not rxn_ids or rxn_ids == '-':
        return pd.DataFrame()
    
    # Split using multiple separators
    pattern = r'[;|,]|' + re.escape(cfg.SPLITER)
    rxn_id_array = [rxn_id.strip() for rxn_id in re.split(pattern, rxn_ids) if rxn_id.strip()]
    
    rxn_list = []
    
    for rxn_id in rxn_id_array:
        try:
            file_id = rxn_id.replace(":", "_") if ':' in rxn_id else rxn_id
            file_path = f"{cfg.DIR_RXN_JSON}{file_id}.json"
            
            try:
                item = ftool.read_json_file(file_path)
                if item:
                    rxn_list.append(item)
                else:
                    print(f"Warning: Empty or invalid JSON read for {file_path}")
            except FileNotFoundError:
                print(f"Warning: Reaction file not found {file_path}")
                
        except Exception as e:
            print(f"Warning: Error processing reaction ID {rxn_id} : {e}")
            continue
    
    if not rxn_list:
        print("Warning: No reaction files were successfully read")
        return pd.DataFrame()
    
    try:
        return pd.json_normalize(rxn_list)
    except Exception as e:
        print(f"Warning: Data normalization failed: {e}")
        return pd.DataFrame()


def get_rxn_details_list(rxn_string: str, rxn_bank: pd.DataFrame, spliter: str = cfg.SPLITER) -> List[Reaction]:
    """
    Parse string containing multiple reaction IDs, returns detailed information list for each reaction
    
    Args:
        rxn_string (str): Reaction string
        rxn_bank (pd.DataFrame): Reaction database DataFrame
        spliter (str, optional): Delimiter, defaults to SPLITER from config file
        
    Returns:
        list: List of valid Reaction objects
    """
    if not rxn_string or rxn_string == '-':
        return []
        
    rxn_ids = [rxn_id.strip() for rxn_id in rxn_string.split(spliter) if rxn_id.strip()]
    return [rxn for rxn in (get_rxn_detail(rxn_id, rxn_bank) for rxn_id in rxn_ids) if rxn is not None]


def get_rxn_details_batch(df_rxns: pd.DataFrame, rxn_bank: pd.DataFrame, 
                         rxn_id_column: str = 'RXNRECer', spliter: str = cfg.SPLITER) -> pd.DataFrame:
    """
    Batch process reaction data in DataFrame, add reaction details for each row
    """
    result_df = df_rxns.copy()
    result_df['RXN_details'] = result_df[rxn_id_column].apply(
        lambda x: get_rxn_details_list(x, rxn_bank, spliter)
    )
    return result_df


def merge_reaction_with_s3_info(RXN_details: List[Reaction], s3_info: List[Dict]) -> List[Dict]:
    """
    Supplement S3 information back to each reaction, generate JSON data for frontend parsing
    """
    if not RXN_details or not s3_info:
        return []
    
    s3_lookup = {item['reaction_id']: item for item in s3_info if 'reaction_id' in item}
    
    merged_reactions = []
    for reaction in RXN_details:
        if reaction is None:
            continue
            
        reaction_dict = reaction.to_dict()
        s3_data = s3_lookup.get(reaction.reaction_id, {})
        
        enriched_reaction = {
            **reaction_dict,
            's3_selected': s3_data.get('selected', 'no'),
            's3_rank': s3_data.get('rank', None),
            's3_confidence': s3_data.get('confidence', 0.0),
            's3_reason': s3_data.get('reason', ''),
            'is_selected': s3_data.get('selected', 'no') == 'yes',
            'selection_rank': s3_data.get('rank', None),
            'confidence_score': s3_data.get('confidence', 0.0),
            'selection_reason': s3_data.get('reason', '')
        }
        merged_reactions.append(enriched_reaction)
    
    return merged_reactions


def create_frontend_friendly_json(RXN_details: List[Reaction], s3_info: List[Dict], output_file: Optional[str] = None) -> Optional[Dict]:
    """
    Create frontend-friendly JSON file containing complete reaction information and S3 scoring
    """
    merged_data = merge_reaction_with_s3_info(RXN_details, s3_info)
    
    frontend_data = {
        'reactions': merged_data,
        'summary': {
            'total_reactions': len(merged_data),
            'selected_reactions': sum(1 for r in merged_data if r.get('is_selected', False)),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'data_source': 'RXNRECer S3 Analysis'
        }
    }
    
    if output_file:
        ftool.write_json_file(frontend_data, output_file)
        print(f"✅ Frontend-friendly JSON file generated: {output_file}")
        return None
    else:
        return frontend_data


def format_obj(x, ndigits: int = 6):
    """Recursively process cell content, preserve floats to specified decimal places"""
    if isinstance(x, (np.floating, float)):
        return round(float(x), ndigits)
    elif isinstance(x, dict):
        return {k: format_obj(v, ndigits) for k, v in x.items()}
    elif isinstance(x, (list, tuple)):
        return [format_obj(v, ndigits) for v in x]
    else:
        return x
    
    
def simplify_rxn_details_fields(rxn_details_list: List[Dict]) -> Tuple[List[str], List[str], List[str]]:
    """
    Extract key fields from reaction details list.
    """
    reaction_ec = []
    reaction_equation = []
    reaction_equation_ref_chebi = []
    
    for item in rxn_details_list:
        rxn_id = item.get('reaction_id', '-')
        
        if rxn_id == '-':
            reaction_ec.append('-')
            reaction_equation.append('-')
            reaction_equation_ref_chebi.append('-')
        else:
            details = item.get('reaction_details', {})
            reaction_ec.append(details.get('reaction_ec', '-'))
            reaction_equation.append(details.get('reaction_equation', '-'))
            reaction_equation_ref_chebi.append(details.get('reaction_equation_ref_chebi', '-'))
            
    return reaction_ec, reaction_equation, reaction_equation_ref_chebi

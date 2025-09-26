"""
Biological utility functions for RXNRECer
"""

import sys
import os
from rxnrecer.config import config as cfg
from rxnrecer.utils import file_utils as ftool
import re
import subprocess
import tempfile
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union
from pandarallel import pandarallel
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


# Removed duplicated dataframe_to_fasta to avoid divergence; use file_utils.dataframe_to_fasta


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

def get_rxn_detail(rxn_id, rxn_bank):
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

    Note:
        Using Reaction object provides unified interface for subsequent processing
    """
    from rxnrecer.lib.rxn.Reaction import Reaction
    
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

def get_rxn_details_from_rxn_json(rxn_ids):
    """
    Get reaction details from JSON files with support for multiple separators and exception handling

    Args:
        rxn_ids (str): Reaction ID string supporting multiple separators (; | ,)
            e.g., 'RHEA:14709;RHEA:24076;RHEA:32187' or 'RHEA:14709|RHEA:24076'

    Returns:
        pd.DataFrame: DataFrame containing reaction details, returns empty DataFrame if error occurs

    Note:
        Automatically handles multiple separators including ; | , etc.
        Exception handling for missing JSON files, skips invalid files
    """
    import pandas as pd
    
    if not rxn_ids or rxn_ids == '-':
        return pd.DataFrame()
    
    # Handle multiple separators including ; | , etc.
    separators = [';', '|', ',', cfg.SPLITER]
    rxn_id_array = []
    
    # Try different separators for splitting
    for sep in separators:
        if sep in rxn_ids:
            rxn_id_array = [rxn_id.strip() for rxn_id in rxn_ids.split(sep) if rxn_id.strip()]
            break
    
    # If no separator found, treat as single ID
    if not rxn_id_array:
        rxn_id_array = [rxn_ids.strip()]
    
    rxn_list = []  # Store JSON data for each reaction
    
    for rxn_id in rxn_id_array:
        try:
            # Handle RHEA: format ID, convert to filename format
            if ':' in rxn_id:
                file_id = rxn_id.replace(":", "_")
            else:
                file_id = rxn_id
            
            # Build file path
            file_path = f"{cfg.DIR_RXN_JSON}{file_id}.json"
            
            # Check if file exists
            if not os.path.exists(file_path):
                print(f"Warning: Reaction file not found {file_path}")
                continue
            
            # Read JSON file
            item = ftool.read_json_file(file_path)
            if item:  # Ensure successful reading
                rxn_list.append(item)
            else:
                print(f"Warning: Unable to read reaction file {file_path}")
                
        except Exception as e:
            print(f"Warning: Error processing reaction ID {rxn_id} : {e}")
            continue
    
    # If no files were successfully read, return empty DataFrame
    if not rxn_list:
        print(f"Warning: No reaction files were successfully read")
        return pd.DataFrame()
    
    try:
        # Use pandas json_normalize to process data
        res = pd.json_normalize(rxn_list)
        return res
    except Exception as e:
        print(f"Warning: Data normalization failed: {e}")
        return pd.DataFrame()


def get_rxn_details_list(rxn_string, rxn_bank, spliter=cfg.SPLITER):
    """
    Parse string containing multiple reaction IDs, returns detailed information list for each reaction
    
    Args:
        rxn_string (str): Reaction string, may contain multiple reaction IDs separated by delimiter
            e.g., 'RHEA:12345|RHEA:67890' or 'RHEA:12345'
        rxn_bank (pd.DataFrame): Reaction database DataFrame
        spliter (str, optional): Delimiter, defaults to SPLITER from config file
        
    Returns:
        list: List of Reaction objects, each object contains complete reaction information
        
    Examples:
        >>> get_rxn_details_list('RHEA:12345|RHEA:67890', rxn_bank, '|')
        [<Reaction object>, <Reaction object>]
        
        >>> get_rxn_details_list('-', rxn_bank)
        []
    """
    # Handle empty values or no reaction cases
    if not rxn_string or rxn_string == '-':
        return []
    else:
        # Split reaction ID string by delimiter
        rxn_ids = [rxn_id.strip() for rxn_id in rxn_string.split(spliter) if rxn_id.strip()]
        # Get detailed information for each reaction ID, filter out None values
        RXN_details = [get_rxn_detail(rxn_id, rxn_bank) for rxn_id in rxn_ids]
        # Filter out None values, only return valid Reaction objects
        return [rxn for rxn in RXN_details if rxn is not None] 


def get_rxn_details_batch(df_rxns, rxn_bank, rxn_id_column='RXNRECer', spliter=cfg.SPLITER):
    """
    Batch process reaction data in DataFrame, add reaction details for each row
    
    Args:
        df_rxns (pd.DataFrame): DataFrame containing reaction data
        rxn_bank (pd.DataFrame): Reaction database DataFrame
        rxn_id_column (str, optional): Reaction ID column name, defaults to'RXNRECer'
        spliter (str, optional): Delimiter, defaults to SPLITER from config file
        
    Returns:
        pd.DataFrame: Expanded DataFrame with new column:
            - RXN_details: List containing all Reaction objects for that row
            
    Examples:
        >>> df = pd.DataFrame({'RXNRECer': ['RHEA:12345', 'RHEA:67890|RHEA:11111']})
        >>> result = get_rxn_details_batch(df, rxn_bank)
        >>> result['RXN_details'][0]  # Reaction details for first row
        [<Reaction object>]
        >>> result['RXN_details'][1]  # Reaction details for second row (contains 2 reactions)
        [<Reaction object>, <Reaction object>]
        
    Note:
        此函数会为每行创建一个RXN_details列，List containing all Reaction objects for that row
        If a row contains multiple reaction IDs (separated by delimiter), returns list containing multiple Reaction objects
        Invalid reaction IDs will be filtered out and will not appear in results
    """
    
    # Create a copy to avoid modifying original
    result_df = df_rxns.copy()
    pandarallel.initialize()
    result_df['RXN_details'] = result_df[rxn_id_column].parallel_apply(
        lambda x: get_rxn_details_list(x, rxn_bank, spliter)
    )
    
    return result_df


def merge_reaction_with_s3_info(RXN_details, s3_info):
    """
    Supplement S3 information back to each reaction, generate JSON data for frontend parsing
    
    Args:
        RXN_details (list): Reaction details list, each element is a Reaction object
        s3_info (list): S3 information list, each element contains reaction_id, selected, rank, confidence, reason
        
    Returns:
        list: Merged reaction information list, each element is a dictionary containing S3 information
        
    Examples:
        >>> RXN_details = [reaction_obj1, reaction_obj2, reaction_obj3]
        >>> s3_info = [
        ...     {'reaction_id': 'RHEA:14709', 'selected': 'yes', 'rank': 1, 'confidence': 0.95, 'reason': '...'},
        ...     {'reaction_id': 'RHEA:24076', 'selected': 'no', 'confidence': 0.2, 'reason': '...'},
        ...     {'reaction_id': 'RHEA:32187', 'selected': 'no', 'confidence': 0.1, 'reason': '...'}
        ... ]
        >>> merged = merge_reaction_with_s3_info(RXN_details, s3_info)
        >>> # Result contains complete reaction information and S3 scoring information
    """
    if not RXN_details or not s3_info:
        return []
    
    # Create S3 information lookup dictionary with reaction_id as key
    s3_lookup = {}
    for s3_item in s3_info:
        if 'reaction_id' in s3_item:
            s3_lookup[s3_item['reaction_id']] = s3_item
    
    merged_reactions = []
    
    for reaction in RXN_details:
        if reaction is None:
            continue
            
        # Get basic reaction information
        reaction_dict = reaction.to_dict()
        
        # Find corresponding S3 information
        s3_data = s3_lookup.get(reaction.reaction_id, {})
        
        # Merge information
        enriched_reaction = {
            # Basic reaction information
            **reaction_dict,
            
            # S3 scoring information
            's3_selected': s3_data.get('selected', 'no'),
            's3_rank': s3_data.get('rank', None),
            's3_confidence': s3_data.get('confidence', 0.0),
            's3_reason': s3_data.get('reason', ''),
            
            # Frontend-friendly fields
            'is_selected': s3_data.get('selected', 'no') == 'yes',
            'selection_rank': s3_data.get('rank', None),
            'confidence_score': s3_data.get('confidence', 0.0),
            'selection_reason': s3_data.get('reason', '')
        }
        
        merged_reactions.append(enriched_reaction)
    
    return merged_reactions


def create_frontend_friendly_json(RXN_details, s3_info, output_file=None):
    """
    Create frontend-friendly JSON file containing complete reaction information and S3 scoring
    
    Args:
        RXN_details (list): Reaction details list
        s3_info (list): S3 information list
        output_file (str, optional): Output file path, returns dictionary if None
        
    Returns:
        dict or None: Returns dictionary if output_file is None, otherwise returns None (writes to file)
        
    Examples:
        >>> # Generate JSON file
        >>> create_frontend_friendly_json(RXN_details, s3_info, 'output.json')
        
        >>> # Return dictionary data
        >>> data = create_frontend_friendly_json(RXN_details, s3_info)
        >>> print(json.dumps(data, indent=2))
    """
    import time
    
    merged_data = merge_reaction_with_s3_info(RXN_details, s3_info)
    
    # Create frontend-friendly data structure
    frontend_data = {
        'reactions': merged_data,
        'summary': {
            'total_reactions': len(merged_data),
            'selected_reactions': len([r for r in merged_data if r.get('is_selected', False)]),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'data_source': 'RXNRECer S3 Analysis'
        }
    }
    
    if output_file:
        # Write to file
        ftool.write_json_file(frontend_data, output_file)
        print(f"✅ Frontend-friendly JSON file generated: {output_file}")
        return None
    else:
        # Return dictionary data
        return frontend_data

def format_obj(x, ndigits=6):
    """递归处理单元格内容，保留浮点数到指定小数位"""
    if isinstance(x, (np.floating, float)):
        return round(float(x), ndigits)
    elif isinstance(x, dict):
        return {k: format_obj(v, ndigits) for k, v in x.items()}
    elif isinstance(x, (list, tuple)):
        return [format_obj(v, ndigits) for v in x]
    else:
        return x
    
    
def simplify_rxn_details_fields(rxn_details_list):
    
    reaction_ec =[]
    reaction_equation = []
    reaction_equation_ref_chebi = []
    
    for item in rxn_details_list:
        rxn_id = item['reaction_id']
        
        if rxn_id == '-':
            rxn_ec = '-'
            rxn_equ = '-'
            rxn_equ_ref_chebi = '-'
        else:
            rxn_ec = item['reaction_details']['reaction_ec']
            rxn_equ = item['reaction_details']['reaction_equation']
            rxn_equ_ref_chebi = item['reaction_details']['reaction_equation_ref_chebi']
            
        reaction_ec.append(rxn_ec)
        reaction_equation.append(rxn_equ)
        reaction_equation_ref_chebi.append(rxn_equ_ref_chebi)
    
    return reaction_ec, reaction_equation, reaction_equation_ref_chebi
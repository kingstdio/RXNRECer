"""
File utility functions for RXNRECer
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import Union, Dict, List, Any
from Bio import SeqIO


def ensure_dir(directory: str) -> None:
    """Ensure directory exists, create if not."""
    Path(directory).mkdir(parents=True, exist_ok=True)


def fasta_to_dataframe(fasta_file: str) -> pd.DataFrame:
    """
    Read FASTA file and convert to DataFrame.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        DataFrame with columns ['uniprot_id', 'seq']
    """
    data = [(record.id, str(record.seq)) for record in SeqIO.parse(fasta_file, "fasta")]
    return pd.DataFrame(data, columns=["uniprot_id", "seq"])


def dataframe_to_fasta(df: pd.DataFrame, output_file: str, 
                      id_col: str = "uniprot_id", seq_col: str = "seq") -> None:
    """
    Convert DataFrame to FASTA file.
    
    Args:
        df: DataFrame with sequence data
        output_file: Output FASTA file path
        id_col: Column name for sequence IDs
        seq_col: Column name for sequences
    """
    with open(output_file, 'w') as f:
        for _, row in df.iterrows():
            f.write(f">{row[id_col]}\n")
            f.write(f"{row[seq_col]}\n")


def read_json_file(file_path: str) -> Dict[str, Any]:
    """
    Read and parse JSON file with data cleaning.
    
    Args:
        file_path: Path to JSON file
        
    Returns:
        Parsed JSON data
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    return clean_data(data)


def write_json_file(data: Dict[str, Any], file_path: str, indent: int = 2) -> None:
    """
    Write data to JSON file.
    
    Args:
        data: Data to write
        file_path: Output file path
        indent: JSON indentation
    """
    ensure_dir(os.path.dirname(file_path))
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=indent)


def clean_string(s: str) -> str:
    """Clean string by removing extra escape characters."""
    if isinstance(s, str):
        return s.replace('\\\\', '\\')
    return s


def clean_data(data: Any) -> Any:
    """
    Recursively clean data structure by removing escape characters.
    
    Args:
        data: Data to clean (dict, list, or primitive type)
        
    Returns:
        Cleaned data
    """
    if isinstance(data, dict):
        return {key: clean_data(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [clean_data(item) for item in data]
    else:
        return clean_string(data)


def save_dataframe(df: pd.DataFrame, output_file: str, 
                  output_format: str = 'tsv') -> None:
    """
    Save DataFrame to file in specified format.
    
    Args:
        df: DataFrame to save
        output_file: Output file path
        output_format: Output format ('tsv', 'csv', 'json', 'feather')
    """
    ensure_dir(os.path.dirname(output_file))
    
    if output_format == 'tsv':
        df.to_csv(output_file, sep='\t', index=False)
    elif output_format == 'csv':
        df.to_csv(output_file, index=False)
    elif output_format == 'json':
        df.to_json(output_file, orient='records', indent=2)
    elif output_format == 'feather':
        df.to_feather(output_file)
    else:
        raise ValueError(f"Unsupported output format: {output_format}")


def load_dataframe(file_path: str, file_format: str = None) -> pd.DataFrame:
    """
    Load DataFrame from file.
    
    Args:
        file_path: Path to file
        file_format: File format (auto-detected if None)
        
    Returns:
        Loaded DataFrame
    """
    if file_format is None:
        file_format = Path(file_path).suffix.lower()
    
    if file_format in ['.tsv', '.txt']:
        return pd.read_csv(file_path, sep='\t')
    elif file_format == '.csv':
        return pd.read_csv(file_path)
    elif file_format == '.json':
        return pd.read_json(file_path)
    elif file_format == '.feather':
        return pd.read_feather(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_format}")


def get_file_hash(file_path: str) -> str:
    """
    Calculate MD5 hash of file content.
    
    Args:
        file_path: Path to file
        
    Returns:
        MD5 hash string
    """
    import hashlib
    
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_dataframe_hash(df: pd.DataFrame) -> str:
    """
    Calculate MD5 hash of DataFrame content.
    
    Args:
        df: DataFrame to hash
        
    Returns:
        MD5 hash string
    """
    import hashlib
    
    df_str = df.to_string(index=False, header=True)
    return hashlib.md5(df_str.encode('utf-8')).hexdigest()

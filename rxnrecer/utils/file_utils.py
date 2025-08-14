"""
File utility functions for RXNRECer
"""

import os
import json
import shutil
import subprocess
import time
import pandas as pd
from pathlib import Path
from typing import  Dict,  Any
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

#region: download file
def downlod(
    download_url: str,
    save_file: str,
    verbos: bool = False,
    connections: int = 16,
    timeout: int = 300,
    allow_overwrite: bool = False,
) -> None:
    """
    Download a file using aria2c if available.

    - Uses multiple connections for speed by default.
    - If aria2c is missing, prints guidance and raises an error.

    Args:
        download_url: Source URL
        save_file: Destination file path
        verbos: Print verbose aria2c output if True
        connections: Parallel connections for aria2c (-x/-s)
        timeout: Per-connection timeout (seconds)
        allow_overwrite: Overwrite existing file if True
    """
    dest_dir = os.path.dirname(save_file) or "."
    ensure_dir(dest_dir)

    aria2c_path = shutil.which('aria2c')
    if not aria2c_path:
        msg = (
            "aria2c not found. Please install aria2 (e.g., `sudo apt install aria2` or "
            "`conda install -c conda-forge aria2`)."
        )
        print(msg)
        raise RuntimeError(msg)

    # If not allowed to overwrite and file already exists, skip
    if not allow_overwrite and os.path.exists(save_file):
        if verbos:
            print(f"Skip existing file: {save_file}")
        return

    file_name = os.path.basename(save_file)
    cmd = [
        aria2c_path,
        f"--allow-overwrite={'true' if allow_overwrite else 'false'}",
        '--auto-file-renaming=false',
        '-x', str(connections),
        '-s', str(connections),
        '-k', '1M',
        '--timeout', str(timeout),
        '-m', '3',
        '-d', dest_dir,
        '-o', file_name,
        download_url,
    ]
    if not verbos:
        cmd.extend(['--console-log-level=warn', '--summary-interval=0'])

    if verbos:
        print(' '.join(cmd))

    subprocess.run(cmd, check=True)


#endregion


#region DataFrame表格转fasta文件
def table2fasta(table, file_out):
    """DataFrame表格转fasta文件, 输入两列，【序列名称，序列】

    Args:
        table (DataFrame): 包含序列名称、序列的DataFame
        file_out (_type_): 输出fasta文件路径
    """
    file = open(file_out, 'w')
    for index, row in table.iterrows():
        file.write(f'>{row.iloc[0]}\n')
        file.write(f'{ row.iloc[1]}\n')
    file.close()
    # print('Write finished')
#endregion
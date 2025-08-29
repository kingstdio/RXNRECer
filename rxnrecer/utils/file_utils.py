"""
File utility functions for RXNRECer
"""
import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/')
import hashlib
import json
import shutil
import subprocess
import pandas as pd
from pathlib import Path
from typing import  Dict,  Any
from Bio import SeqIO
import hashlib
from rxnrecer.config import config as cfg

def get_project_root() -> Path:
    """
    获取项目根目录，优先使用环境变量，否则智能查找
    
    Returns:
        Path: 项目根目录路径
        
    Note:
        优先检查环境变量RXNRECER_PROJECT_ROOT，如果没有设置，
        则智能查找包含data/和ckpt/目录的父目录
    """
    # 优先使用环境变量
    env_root = os.environ.get('RXNRECER_PROJECT_ROOT')
    if env_root:
        env_path = Path(env_root)
        if (env_path / 'data').exists() and (env_path / 'ckpt').exists():
            return env_path
        else:
            print(f"⚠️  警告: 环境变量RXNRECER_PROJECT_ROOT指定的路径 {env_root} 不包含必要的data和ckpt目录")
    
    # 智能查找项目根目录
    # 从当前工作目录开始，向上查找包含data和ckpt的目录
    current_dir = Path.cwd()
    
    # 检查当前目录
    if (current_dir / 'data').exists() and (current_dir / 'ckpt').exists():
        return current_dir
    
    # 向上查找父目录
    for parent in current_dir.parents:
        if (parent / 'data').exists() and (parent / 'ckpt').exists():
            return parent
    
    # 如果都找不到，尝试查找脚本文件所在目录
    script_dir = Path(__file__).parent.parent.parent
    if (script_dir / 'data').exists() and (script_dir / 'ckpt').exists():
        return script_dir
    
    # 最后尝试脚本目录的父目录
    script_parent = script_dir.parent
    if (script_parent / 'data').exists() and (script_parent / 'ckpt').exists():
        return script_parent
    
    # 如果都找不到，使用当前目录并给出警告
    print(f"⚠️  警告: 无法找到包含data和ckpt目录的项目根目录")
    print(f"   当前工作目录: {current_dir}")
    print(f"   脚本目录: {script_dir}")
    print(f"   请确保在正确的项目目录下运行，或设置RXNRECER_PROJECT_ROOT环境变量")
    return current_dir


def get_data_root() -> str:
    """获取数据目录路径"""
    return str(get_project_root() / 'data')


def get_ckpt_root() -> str:
    """获取模型检查点目录路径"""
    return str(get_project_root() / 'ckpt')


def get_results_root() -> str:
    """获取结果输出目录路径"""
    return str(get_project_root() / 'results')


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


def write_json_file(data: Any, file_path: str, indent: int = 2, overwrite: bool = True) -> None:
    """
    Write data to JSON file with automatic directory creation.
    
    Args:
        data: Data to write (dict, list, or any JSON-serializable object)
        file_path: Output file path
        indent: JSON indentation
        overwrite: Whether to overwrite existing file (default: True)
    """
    # Ensure directory exists
    ensure_dir(os.path.dirname(file_path))
    
    # Check if file exists and handle accordingly
    if os.path.exists(file_path) and not overwrite:
        print(f'File exists and overwrite=False: {file_path}')
        return
    
    # Write data to file
    with open(file_path, 'w') as f:
        if hasattr(data, 'to_json'):
            # Handle objects with to_json method (like Reaction objects)
            f.write(data.to_json())
        else:
            # Handle regular data types
            json.dump(data, f, indent=indent, ensure_ascii=False)


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


def read_json_as_object(file_path: str) -> object:
    """
    Read JSON file and convert to object with dot notation support.
    
    Args:
        file_path: Path to JSON file
        
    Returns:
        Object with dot notation access to JSON data
        
    Example:
        prompt_obj = read_json_as_object('prompts.json')
        text = prompt_obj.prompts1.text
    """
    from types import SimpleNamespace
    
    data = read_json_file(file_path)
    return json.loads(json.dumps(data), object_hook=lambda d: SimpleNamespace(**d))



def get_cache_filename(input_file, mode, output_format):
    """生成缓存文件名，基于输入文件内容、模式和输出格式的hash"""
    try:
        with open(input_file, 'rb') as f:
            file_content = f.read()
        # 使用文件内容、模式和输出格式生成hash
        content_hash = hashlib.md5(f'{file_content}_{mode}_{output_format}'.encode()).hexdigest()
        return f"cache_{content_hash}"
    except Exception as e:
        print(f"Warning: Could not read file for cache key: {e}")
        return None


def check_cache(cache_filename):
    """检查缓存是否存在，如果存在则返回缓存文件路径"""
    
    cache_file = f'{cfg.CACHE_DIR}{cache_filename}.pkl'
    # 创建缓存目录
    
    # 检查缓存文件是否存在
    if os.path.exists(cache_file):
        print(f"📋 Found cached results: {cache_filename}")
        return True
    else:
        if not os.path.exists(f'{cfg.CACHE_DIR}'):
            os.makedirs(f'{cfg.CACHE_DIR}', exist_ok=True)
        return False
    


def save_to_cache(cache_data,cache_filename):
    """保存结果到缓存"""
    try:
        cache_file = f'{cfg.CACHE_DIR}{cache_filename}.pkl'
        cache_data.to_pickle(cache_file)
    except Exception as e:
        print(f"Warning: Failed to cache results: {e}")
        
def load_from_cache(cache_filename):
    """从缓存加载结果"""
    try:
        cache_file = f'{cfg.CACHE_DIR}{cache_filename}.pkl'
        return pd.read_pickle(cache_file)
    except Exception as e:
        print(f"Warning: Failed to load cached results: {e}")
        return None


def clear_cache(cache_filename=None, older_than_days=None):
    """
    清理缓存文件
    
    Args:
        cache_filename (str, optional): 指定要删除的缓存文件
        older_than_days (int, optional): 删除指定天数前的缓存文件
    
    Returns:
        int: 删除的文件数量
    """
    try:
        deleted_count = 0
        
        if cache_filename:
            # 删除指定文件
            cache_file = f'{cfg.CACHE_DIR}{cache_filename}.pkl'
            if os.path.exists(cache_file):
                os.remove(cache_file)
                deleted_count += 1
                print(f"🗑️  Deleted cache file: {cache_filename}")
            else:
                print(f"Warning: Cache file not found: {cache_filename}")
        elif older_than_days:
            # 删除指定天数前的文件
            import time
            current_time = time.time()
            cutoff_time = current_time - (older_than_days * 24 * 3600)
            
            for filename in os.listdir(cfg.CACHE_DIR):
                if filename.endswith('.pkl'):
                    file_path = os.path.join(cfg.CACHE_DIR, filename)
                    if os.path.getmtime(file_path) < cutoff_time:
                        os.remove(file_path)
                        deleted_count += 1
                        print(f"🗑️  Deleted old cache file: {filename}")
        else:
            # 删除所有缓存文件
            for filename in os.listdir(cfg.CACHE_DIR):
                if filename.endswith('.pkl'):
                    file_path = os.path.join(cfg.CACHE_DIR, filename)
                    os.remove(file_path)
                    deleted_count += 1
            
            if deleted_count > 0:
                print(f"🗑️  Deleted {deleted_count} cache files")
        
        return deleted_count
        
    except PermissionError:
        print(f"Error: Permission denied clearing cache")
        return 0
    except OSError as e:
        print(f"Error: OS error clearing cache: {e}")
        return 0
    except Exception as e:
        print(f"Warning: Unexpected error clearing cache: {e}")
        return 0


def get_cache_info():
    """
    获取缓存信息
    
    Returns:
        dict: 包含缓存统计信息的字典
    """
    try:
        if not os.path.exists(cfg.CACHE_DIR):
            return {
                'cache_dir': cfg.CACHE_DIR,
                'total_files': 0,
                'total_size_mb': 0.0,
                'oldest_file': None,
                'newest_file': None
            }
        
        cache_files = [f for f in os.listdir(cfg.CACHE_DIR) if f.endswith('.pkl')]
        
        if not cache_files:
            return {
                'cache_dir': cfg.CACHE_DIR,
                'total_files': 0,
                'total_size_mb': 0.0,
                'oldest_file': None,
                'newest_file': None
            }
        
        total_size = 0
        file_times = []
        
        for filename in cache_files:
            file_path = os.path.join(cfg.CACHE_DIR, filename)
            file_size = os.path.getsize(file_path)
            file_time = os.path.getmtime(file_path)
            
            total_size += file_size
            file_times.append((filename, file_time))
        
        # 按时间排序
        file_times.sort(key=lambda x: x[1])
        
        return {
            'cache_dir': cfg.CACHE_DIR,
            'total_files': len(cache_files),
            'total_size_mb': round(total_size / (1024 * 1024), 2),
            'oldest_file': file_times[0][0] if file_times else None,
            'newest_file': file_times[-1][0] if file_times else None
        }
        
    except Exception as e:
        print(f"Warning: Failed to get cache info: {e}")
        return None
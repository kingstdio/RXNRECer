"""
RXNRECer: Deep learning framework for predicting enzyme-catalyzed reactions from protein sequences.

This package provides a comprehensive framework for enzyme reaction prediction using
state-of-the-art protein language models and deep learning techniques.

Author: Zhenkun Shi
Email: zhenkun.shi@tib.cas.cn
Project: https://github.com/kingstdio/RXNRECer
"""

import os
import hashlib
import json
from pathlib import Path

# Version information
__version__ = "1.3.1"
__version_info__ = (1, 3, 1)

def get_version():
    """Get the version string."""
    return __version__

def get_version_info():
    """Get the version tuple."""
    return __version_info__

def get_full_version():
    """Get the full version string with build info."""
    return f"{__version__} (2024-09-26)"

__author__ = "Zhenkun Shi"
__email__ = "zhenkun.shi@tib.cas.cn"
__project__ = "RXNRECer"
__url__ = "https://github.com/kingstdio/RXNRECer"

# Import main components
from .config import config
from .cli.predict import main as predict_main

def check_data_files():
    """
    Check if required data files are available.
    
    Returns:
        bool: True if all required files exist, False otherwise
    """
    required_dirs = [
        "data/sample",
        "ckpt/rxnrecer",
        "ckpt/prostt5",
        "extools"
    ]
    
    missing_dirs = []
    for dir_path in required_dirs:
        if not os.path.exists(dir_path) or not os.listdir(dir_path):
            missing_dirs.append(dir_path)
    
    if missing_dirs:
        print("⚠️  Required data files are missing!")
        print("📁 Missing directories:")
        for dir_path in missing_dirs:
            print(f"   - {dir_path}")
        print("\n📥 To download data files, run:")
        print("   rxnrecer-download-data")
        print("\n📚 Or manually download from:")
        print("   - Data (国内镜像): http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/data.tar.gz")
        print("   - Data (AWS镜像): https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz")
        print("   - Models (国内镜像): http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/ckpt.tar.gz")
        print("   - Models (AWS镜像): https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz")
        print("   - Extools (国内镜像): http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/extools.tar.gz")
        print("   - Extools (AWS镜像): https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/extools.tar.gz")
        return False
    
    print("✅ All required data files are available")
    return True

def select_fastest_mirror(urls, timeout=5):
    """
    Select the fastest available mirror by testing connection speed.
    
    Args:
        urls (list): List of mirror URLs to test
        timeout (int): Connection timeout in seconds
    
    Returns:
        str: Fastest available URL, or first URL if all fail
    """
    import socket
    import time
    from urllib.parse import urlparse
    
    fastest_url = urls[0]  # Default to first URL
    fastest_time = float('inf')
    
    print("🔍 Testing mirror speeds...")
    
    for url in urls:
        try:
            parsed = urlparse(url)
            host = parsed.hostname
            port = parsed.port or (443 if parsed.scheme == 'https' else 80)
            
            # Test connection speed
            start_time = time.time()
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(timeout)
            result = sock.connect_ex((host, port))
            sock.close()
            
            if result == 0:
                connect_time = time.time() - start_time
                print(f"   {host}: {connect_time:.2f}s")
                
                if connect_time < fastest_time:
                    fastest_time = connect_time
                    fastest_url = url
            else:
                print(f"   {host}: Connection failed")
                
        except Exception as e:
            print(f"   {parsed.hostname}: Error - {e}")
            continue
    
    print(f"✅ Selected fastest mirror: {urlparse(fastest_url).hostname}")
    return fastest_url

def download_data_files(force=False, data_only=False, models_only=False, extools_only=False):
    """
    Download required data files automatically.
    
    Args:
        force (bool): Force download even if files exist
        data_only (bool): Download only data files
        models_only (bool): Download only model files
        extools_only (bool): Download only external tools
    
    Returns:
        bool: True if download successful, False otherwise
    """
    try:
        import subprocess
        import urllib.request
        import tarfile
        
        print("🚀 Starting automatic data download...")
        
        # Check if wget is available
        try:
            subprocess.run(["wget", "--version"], capture_output=True, check=True)
            use_wget = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            use_wget = False
        
        # Download data files - multiple mirrors for better availability
        data_urls = [
            "http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/data.tar.gz",  # 国内镜像
            "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz"  # AWS镜像
        ]
        ckpt_urls = [
            "http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/ckpt.tar.gz",  # 国内镜像
            "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz"  # AWS镜像
        ]
        extools_urls = [
            "http://s3.biodesign.ac.cn/bucket-rxnrecer-web/files/extools.tar.gz",  # 国内镜像
            "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/extools.tar.gz"  # AWS镜像
        ]
        
        if (not data_only and not models_only and not extools_only) or data_only:
            if not os.path.exists("data/sample") or force:
                print("📥 Downloading data files (~8.8GB)...")
                data_url = select_fastest_mirror(data_urls)
                if use_wget:
                    subprocess.run(["wget", "-O", "data.tar.gz", data_url], check=True)
                else:
                    urllib.request.urlretrieve(data_url, "data.tar.gz")
                
                print("📦 Extracting data files...")
                with tarfile.open("data.tar.gz", "r:gz") as tar:
                    tar.extractall(".")
                os.remove("data.tar.gz")
                print("✅ Data files downloaded and extracted")
        
        if (not data_only and not models_only and not extools_only) or models_only:
            if not os.path.exists("ckpt/rxnrecer") or force:
                print("📥 Downloading model files (~14GB)...")
                ckpt_url = select_fastest_mirror(ckpt_urls)
                if use_wget:
                    subprocess.run(["wget", "-O", "ckpt.tar.gz", ckpt_url], check=True)
                else:
                    urllib.request.urlretrieve(ckpt_url, "ckpt.tar.gz")
                
                print("📦 Extracting model files...")
                with tarfile.open("ckpt.tar.gz", "r:gz") as tar:
                    tar.extractall(".")
                os.remove("ckpt.tar.gz")
                print("✅ Model files downloaded and extracted")
        
        if (not data_only and not models_only and not extools_only) or extools_only:
            if not os.path.exists("extools") or force:
                print("📥 Downloading external tools (~13GB)...")
                extools_url = select_fastest_mirror(extools_urls)
                if use_wget:
                    subprocess.run(["wget", "-O", "extools.tar.gz", extools_url], check=True)
                else:
                    urllib.request.urlretrieve(extools_url, "extools.tar.gz")
                
                print("📦 Extracting external tools...")
                with tarfile.open("extools.tar.gz", "r:gz") as tar:
                    tar.extractall(".")
                os.remove("extools.tar.gz")
                print("✅ External tools downloaded and extracted")
        
        print("🎉 All files downloaded successfully!")
        return True
        
    except Exception as e:
        print(f"❌ Download failed: {e}")
        print("📚 Please download manually following the instructions in data/README.md")
        return False

def get_cache_key(input_file, mode, format, batch_size):
    """
    Generate a cache key for the given input parameters.
    
    Args:
        input_file (str): Input FASTA file path
        mode (str): Prediction mode
        format (str): Output format
        batch_size (int): Batch size
    
    Returns:
        str: Cache key
    """
    # Read file content and generate hash
    with open(input_file, 'rb') as f:
        content = f.read()
    
    # Generate hash from content and parameters
    hash_input = f"{content}{mode}{format}{batch_size}".encode()
    return hashlib.md5(hash_input).hexdigest()

def get_cache_path(cache_key, output_format):
    """
    Get the cache file path for the given cache key.
    
    Args:
        cache_key (str): Cache key
        output_format (str): Output format (tsv/json)
    
    Returns:
        str: Cache file path
    """
    cache_dir = Path("results/cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    extension = "json" if output_format == "json" else "tsv"
    return cache_dir / f"{cache_key}.{extension}"

def is_cached(input_file, mode, format, batch_size):
    """
    Check if results are cached for the given parameters.
    
    Args:
        input_file (str): Input FASTA file path
        mode (str): Prediction mode
        format (str): Output format
        batch_size (int): Batch size
    
    Returns:
        bool: True if cached results exist, False otherwise
    """
    try:
        cache_key = get_cache_key(input_file, mode, format, batch_size)
        cache_path = get_cache_path(cache_key, format)
        return cache_path.exists()
    except:
        return False

def get_cached_result(input_file, mode, format, batch_size):
    """
    Get cached result if available.
    
    Args:
        input_file (str): Input FASTA file path
        mode (str): Prediction mode
        format (str): Output format
        batch_size (int): Batch size
    
    Returns:
        str or None: Cached result content if available, None otherwise
    """
    try:
        cache_key = get_cache_key(input_file, mode, format, batch_size)
        cache_path = get_cache_path(cache_key, format)
        
        if cache_path.exists():
            with open(cache_path, 'r') as f:
                return f.read()
        return None
    except:
        return None

def save_to_cache(input_file, mode, format, batch_size, result_content):
    """
    Save result to cache.
    
    Args:
        input_file (str): Input FASTA file path
        mode (str): Prediction mode
        format (str): Output format
        batch_size (int): Batch size
        result_content (str): Result content to cache
    """
    try:
        cache_key = get_cache_key(input_file, mode, format, batch_size)
        cache_path = get_cache_path(cache_key, format)
        
        with open(cache_path, 'w') as f:
            f.write(result_content)
        
        print(f"💾 Results cached for future use: {cache_path}")
    except Exception as e:
        print(f"⚠️  Failed to cache results: {e}")

# Main prediction function with caching
def predict(input_file, output_file, mode="s1", format="tsv", batch_size=100, use_cache=True):
    """
    Main prediction function for RXNRECer with caching support.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output file
        mode (str): Prediction mode ('s1', 's2', or 's3')
        format (str): Output format ('tsv' or 'json')
        batch_size (int): Batch size for processing
        use_cache (bool): Whether to use caching
    
    Returns:
        bool: True if prediction successful, False otherwise
    """
    try:
        # Check data files first
        if not check_data_files():
            return False
        
        # Check cache if enabled
        if use_cache and is_cached(input_file, mode, format, batch_size):
            print("📋 Using cached results...")
            cached_result = get_cached_result(input_file, mode, format, batch_size)
            if cached_result:
                # Write cached result to output file
                with open(output_file, 'w') as f:
                    f.write(cached_result)
                print(f"✅ Results loaded from cache and saved to {output_file}")
                return True
        
        # Run prediction
        print("🚀 Running prediction...")
        predict_main(input_file, output_file, mode, format, batch_size)
        
        # Cache results if enabled
        if use_cache:
            try:
                with open(output_file, 'r') as f:
                    result_content = f.read()
                save_to_cache(input_file, mode, format, batch_size, result_content)
            except:
                pass
        
        return True
    except Exception as e:
        print(f"❌ Prediction failed: {e}")
        return False

# CLI entry point
def cli():
    """Command line interface entry point."""
    from .cli.predict import main
    main()

# Package information
__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__project__",
    "__url__",
    "config",
    "predict",
    "predict_main",
    "cli",
    "check_data_files",
    "download_data_files",
    "is_cached",
    "get_cached_result",
    "save_to_cache"
]

#!/usr/bin/env python3
"""
Cache management command for RXNRECer.
"""

import click
import os
import sys
from pathlib import Path
import shutil

@click.group()
def cache():
    """Manage RXNRECer prediction cache."""
    pass

def main():
    """Main entry point for cache command."""
    cache()

@cache.command()
def status():
    """Show cache status and statistics."""
    try:
        from rxnrecer import is_cached
        cache_dir = Path("results/cache")
        
        if not cache_dir.exists():
            print("📋 Cache directory does not exist")
            return
        
        cache_files = list(cache_dir.glob("*"))
        total_size = sum(f.stat().st_size for f in cache_files if f.is_file())
        
        print("📋 RXNRECer Cache Status")
        print("=" * 30)
        print(f"Cache directory: {cache_dir.absolute()}")
        print(f"Total files: {len(cache_files)}")
        print(f"Total size: {total_size / 1024 / 1024:.2f} MB")
        
        if cache_files:
            print("\n📁 Cached files:")
            for cache_file in sorted(cache_files):
                size = cache_file.stat().st_size / 1024  # KB
                print(f"   - {cache_file.name}: {size:.1f} KB")
        
    except ImportError:
        print("❌ Error: Could not import RXNRECer. Please install it first:")
        print("   pip install rxnrecer")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

@cache.command()
@click.option('--all', is_flag=True, help='Clear all cache files')
@click.option('--file', type=str, help='Clear specific cache file')
@click.option('--older-than', type=int, help='Clear cache files older than N days')
def clear(all, file, older_than):
    """Clear cache files."""
    try:
        cache_dir = Path("results/cache")
        
        if not cache_dir.exists():
            print("📋 Cache directory does not exist")
            return
        
        if file:
            # Clear specific file
            cache_file = cache_dir / file
            if cache_file.exists():
                cache_file.unlink()
                print(f"🗑️  Cleared cache file: {file}")
            else:
                print(f"❌ Cache file not found: {file}")
        elif all:
            # Clear all cache
            shutil.rmtree(cache_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            print("🗑️  Cleared all cache files")
        elif older_than:
            # Clear old cache files
            import time
            current_time = time.time()
            cutoff_time = current_time - (older_than * 24 * 3600)
            
            cleared_count = 0
            for cache_file in cache_dir.glob("*"):
                if cache_file.is_file() and cache_file.stat().st_mtime < cutoff_time:
                    cache_file.unlink()
                    cleared_count += 1
            
            print(f"🗑️  Cleared {cleared_count} cache files older than {older_than} days")
        else:
            print("❌ Please specify what to clear:")
            print("   --all: Clear all cache")
            print("   --file FILENAME: Clear specific file")
            print("   --older-than N: Clear files older than N days")
            
    except ImportError:
        print("❌ Error: Could not import RXNRECer. Please install it first:")
        print("   pip install rxnrecer")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

@cache.command()
def info():
    """Show detailed cache information."""
    try:
        cache_dir = Path("results/cache")
        
        if not cache_dir.exists():
            print("📋 Cache directory does not exist")
            return
        
        print("📋 RXNRECer Cache Information")
        print("=" * 40)
        print(f"Cache directory: {cache_dir.absolute()}")
        print(f"Cache purpose: Store prediction results to avoid recomputation")
        print(f"Cache key: MD5 hash of input file + parameters")
        print(f"Cache location: results/cache/")
        print("\n💡 Benefits:")
        print("   - Faster repeated predictions")
        print("   - Save computational resources")
        print("   - Automatic cache management")
        print("\n🔧 Usage:")
        print("   - Enable: rxnrecer -i input.fasta -o output.tsv (default)")
        print("   - Disable: rxnrecer -i input.fasta -o output.tsv --no-cache")
        print("   - Status: rxnrecer-cache status")
        print("   - Clear: rxnrecer-cache clear --all")
        
    except ImportError:
        print("❌ Error: Could not import RXNRECer. Please install it first:")
        print("   pip install rxnrecer")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    cache()

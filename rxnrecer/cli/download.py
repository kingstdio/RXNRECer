#!/usr/bin/env python3
"""
Data download command for RXNRECer.
"""

import click
import os
import sys
from pathlib import Path

@click.command()
@click.option('--force', is_flag=True, help='Force download even if files exist')
@click.option('--data-only', is_flag=True, help='Download only data files')
@click.option('--models-only', is_flag=True, help='Download only model files')
@click.option('--output-dir', default='.', help='Output directory for downloads')
def download_data(force, data_only, models_only, output_dir):
    """
    Download RXNRECer data and model files.
    
    This command will download the required data files (~8.6GB) and model files (~11.9GB)
    from AWS S3 to your local machine.
    """
    try:
        # Change to output directory
        os.chdir(output_dir)
        
        # Import download function
        from rxnrecer import download_data_files
        
        print("🚀 RXNRECer Data Download Tool")
        print("=" * 40)
        
        if data_only:
            print("📥 Downloading data files only...")
        elif models_only:
            print("📥 Downloading model files only...")
        else:
            print("📥 Downloading all required files...")
        
        # Run download
        success = download_data_files(force=force)
        
        if success:
            print("\n🎉 Download completed successfully!")
            print("\n📁 Files downloaded to:")
            print(f"   - Data: {os.path.abspath('data/')}")
            print(f"   - Models: {os.path.abspath('ckpt/')}")
            print("\n🚀 You can now use RXNRECer:")
            print("   rxnrecer -i input.fasta -o output.tsv -m s1")
        else:
            print("\n❌ Download failed. Please check the error messages above.")
            sys.exit(1)
            
    except ImportError:
        print("❌ Error: Could not import RXNRECer. Please install it first:")
        print("   pip install rxnrecer")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)

if __name__ == '__main__':
    download_data()

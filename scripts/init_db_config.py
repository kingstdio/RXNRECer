#!/usr/bin/env python
"""
Script to initialize database configuration
Usage: python scripts/init_db_config.py
"""
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from rxnrecer.config.init_db_config import main

if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Documentation version updater for RXNRECer.
Updates version numbers in all Markdown documentation files.
"""

import os
import re
from pathlib import Path
from version import __version__

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent

# Markdown files that contain version references
MD_FILES = [
    "README.md",
    "docs/INSTALL.md", 
    "docs/RELEASE_NOTES.md",
    "ckpt/README.md",
    "data/README.md"
]

# Version patterns to update
VERSION_PATTERNS = [
    # Main version references
    (r'\*\*RXNRECer v\d+\.\d+\.\d+\*\*', f'**RXNRECer v{__version__}**'),
    (r'RXNRECer v\d+\.\d+\.\d+', f'RXNRECer v{__version__}'),
    (r'Version \d+\.\d+\.\d+', f'Version {__version__}'),
    (r'v\d+\.\d+\.\d+', f'v{__version__}'),
    
    # Installation commands
    (r'pip install rxnrecer==\d+\.\d+\.\d+', f'pip install rxnrecer=={__version__}'),
    (r'pip install rxnrecer@\d+\.\d+\.\d+', f'pip install rxnrecer@{__version__}'),
    
    # PyPI links
    (r'https://pypi\.org/project/rxnrecer/\d+\.\d+\.\d+/', f'https://pypi.org/project/rxnrecer/{__version__}/'),
]

def update_md_file(file_path):
    """Update version references in a Markdown file"""
    if not file_path.exists():
        print(f"⚠️  File not found: {file_path}")
        return False
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    changes_made = False
    
    # Apply all version patterns
    for pattern, replacement in VERSION_PATTERNS:
        new_content = re.sub(pattern, replacement, content)
        if new_content != content:
            changes_made = True
            content = new_content
    
    if changes_made:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"✅ Updated {file_path}")
        return True
    else:
        print(f"ℹ️  No changes needed in {file_path}")
        return False

def update_all_docs():
    """Update all documentation files"""
    print(f"🔄 Updating documentation to version {__version__}...")
    
    updated_files = []
    for md_file in MD_FILES:
        file_path = PROJECT_ROOT / md_file
        if update_md_file(file_path):
            updated_files.append(md_file)
    
    if updated_files:
        print(f"\n✅ Updated {len(updated_files)} files:")
        for file in updated_files:
            print(f"   - {file}")
    else:
        print("\nℹ️  No files needed updating")
    
    return updated_files

def main():
    """Main function"""
    print(f"📚 RXNRECer Documentation Updater")
    print(f"Current version: {__version__}")
    print("=" * 50)
    
    updated_files = update_all_docs()
    
    if updated_files:
        print(f"\n💡 Don't forget to commit the changes:")
        print(f"   git add {' '.join(updated_files)}")
        print(f"   git commit -m 'docs: Update version references to {__version__}'")

if __name__ == '__main__':
    main()

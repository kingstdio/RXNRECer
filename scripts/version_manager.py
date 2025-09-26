#!/usr/bin/env python3
"""
Version management script for RXNRECer.
This script provides commands to manage version numbers across the project.
"""

import os
import re
import sys
import argparse
from pathlib import Path
from datetime import datetime

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent

# Files that contain version numbers
VERSION_FILES = [
    "version.py",
    "rxnrecer/__init__.py", 
    "README.md",
    "rxnrecer/cli/predict.py",
    "scripts/build_and_release.py",
    "docs/RELEASE_NOTES.md",
    "docs/INSTALL.md",
    "ckpt/README.md",
    "data/README.md"
]

# Markdown-specific version patterns
MD_VERSION_PATTERNS = [
    # Main version references
    (r'\*\*RXNRECer v\d+\.\d+\.\d+\*\*', f'**RXNRECer v{{version}}**'),
    (r'RXNRECer v\d+\.\d+\.\d+', f'RXNRECer v{{version}}'),
    (r'Version \d+\.\d+\.\d+', f'Version {{version}}'),
    (r'v\d+\.\d+\.\d+', f'v{{version}}'),
    
    # Installation commands
    (r'pip install rxnrecer==\d+\.\d+\.\d+', f'pip install rxnrecer=={{version}}'),
    (r'pip install rxnrecer@\d+\.\d+\.\d+', f'pip install rxnrecer@{{version}}'),
    
    # PyPI links
    (r'https://pypi\.org/project/rxnrecer/\d+\.\d+\.\d+/', f'https://pypi.org/project/rxnrecer/{{version}}/'),
]

def get_current_version():
    """Get current version from version.py"""
    version_file = PROJECT_ROOT / "version.py"
    with open(version_file, 'r') as f:
        content = f.read()
        match = re.search(r'__version__ = "([^"]+)"', content)
        if match:
            return match.group(1)
    return None

def update_version_file(new_version):
    """Update version.py with new version"""
    version_file = PROJECT_ROOT / "version.py"
    
    with open(version_file, 'r') as f:
        content = f.read()
    
    # Update version
    content = re.sub(r'__version__ = "[^"]+"', f'__version__ = "{new_version}"', content)
    
    # Update version info tuple
    major, minor, patch = new_version.split('.')
    content = re.sub(
        r'__version_info__ = \(\d+, \d+, \d+\)',
        f'__version_info__ = ({major}, {minor}, {patch})',
        content
    )
    
    # Update build date
    build_date = datetime.now().strftime("%Y-%m-%d")
    content = re.sub(r'BUILD_DATE = "[^"]+"', f'BUILD_DATE = "{build_date}"', content)
    
    with open(version_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Updated {version_file} to version {new_version}")

def update_readme_version(new_version):
    """Update README.md version references"""
    readme_file = PROJECT_ROOT / "README.md"
    
    with open(readme_file, 'r') as f:
        content = f.read()
    
    # Update main version reference
    content = re.sub(
        r'\*\*RXNRECer v\d+\.\d+\.\d+\*\*',
        f'**RXNRECer v{new_version}**',
        content
    )
    
    with open(readme_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Updated README.md to version {new_version}")

def update_predict_version(new_version):
    """Update predict.py version reference"""
    predict_file = PROJECT_ROOT / "rxnrecer/cli/predict.py"
    
    with open(predict_file, 'r') as f:
        content = f.read()
    
    content = re.sub(
        r'print\(f"RXNRECer v\d+\.\d+\.\d+ - Enzyme Reaction Prediction"\)',
        f'print(f"RXNRECer v{new_version} - Enzyme Reaction Prediction")',
        content
    )
    
    with open(predict_file, 'w') as f:
        f.write(content)
    
    print(f"✅ Updated predict.py to version {new_version}")

def update_md_file_version(file_path, new_version):
    """Update version references in a Markdown file"""
    if not file_path.exists():
        print(f"⚠️  File not found: {file_path}")
        return False
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    changes_made = False
    
    # Apply all MD version patterns
    for pattern, replacement_template in MD_VERSION_PATTERNS:
        replacement = replacement_template.format(version=new_version)
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

def update_all_versions(new_version):
    """Update all version references across the project"""
    print(f"🔄 Updating all version references to {new_version}...")
    
    # Update core version file
    update_version_file(new_version)
    
    # Update other files
    update_readme_version(new_version)
    update_predict_version(new_version)
    
    # Update all Markdown files
    md_files = [f for f in VERSION_FILES if f.endswith('.md')]
    for md_file in md_files:
        file_path = PROJECT_ROOT / md_file
        update_md_file_version(file_path, new_version)
    
    print(f"✅ All version references updated to {new_version}")

def validate_version(version):
    """Validate version format (semantic versioning)"""
    pattern = r'^\d+\.\d+\.\d+$'
    if not re.match(pattern, version):
        raise ValueError(f"Invalid version format: {version}. Expected format: X.Y.Z")
    return True

def bump_version(part):
    """Bump version by part (major, minor, patch)"""
    current = get_current_version()
    if not current:
        raise ValueError("Could not determine current version")
    
    major, minor, patch = map(int, current.split('.'))
    
    if part == 'major':
        major += 1
        minor = 0
        patch = 0
    elif part == 'minor':
        minor += 1
        patch = 0
    elif part == 'patch':
        patch += 1
    else:
        raise ValueError(f"Invalid version part: {part}. Use major, minor, or patch")
    
    new_version = f"{major}.{minor}.{patch}"
    return new_version

def main():
    parser = argparse.ArgumentParser(description="RXNRECer Version Manager")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Show current version
    subparsers.add_parser('current', help='Show current version')
    
    # Set specific version
    set_parser = subparsers.add_parser('set', help='Set specific version')
    set_parser.add_argument('version', help='Version to set (e.g., 1.4.0)')
    
    # Bump version
    bump_parser = subparsers.add_parser('bump', help='Bump version')
    bump_parser.add_argument('part', choices=['major', 'minor', 'patch'], 
                           help='Version part to bump')
    
    args = parser.parse_args()
    
    if args.command == 'current':
        version = get_current_version()
        if version:
            print(f"Current version: {version}")
        else:
            print("❌ Could not determine current version")
            sys.exit(1)
    
    elif args.command == 'set':
        try:
            validate_version(args.version)
            update_all_versions(args.version)
            print(f"\n🎉 Version updated to {args.version}")
            print("💡 Don't forget to commit and tag the changes:")
            print(f"   git add .")
            print(f"   git commit -m 'Bump version to {args.version}'")
            print(f"   git tag v{args.version}")
        except ValueError as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
    
    elif args.command == 'bump':
        try:
            new_version = bump_version(args.part)
            update_all_versions(new_version)
            print(f"\n🎉 Version bumped to {new_version}")
            print("💡 Don't forget to commit and tag the changes:")
            print(f"   git add .")
            print(f"   git commit -m 'Bump version to {new_version}'")
            print(f"   git tag v{new_version}")
        except ValueError as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
    
    else:
        parser.print_help()

if __name__ == '__main__':
    main()

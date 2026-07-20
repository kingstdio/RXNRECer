#!/usr/bin/env python3
"""
Version synchronization script for RXNRECer.
Ensures all version references are consistent across the project.
"""

import sys
import subprocess
from pathlib import Path

# Add project root to Python path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from version import __version__

# Remove duplicate PROJECT_ROOT definition

PYTHON = sys.executable


def run_command(cmd, capture_output=True):
    """Run a command and return output"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=capture_output, text=True)
        if result.returncode != 0:
            print(f"❌ Command failed: {cmd}")
            print(f"Error: {result.stderr}")
            return None
        return result.stdout.strip() if capture_output else True
    except Exception as e:
        print(f"❌ Error running command: {e}")
        return None

def check_git_status():
    """Check if we're in a git repository and get status"""
    if not (PROJECT_ROOT / ".git").exists():
        print("❌ Not in a git repository")
        return False
    
    # Check if there are uncommitted changes
    status = run_command("git status --porcelain")
    if status:
        print("⚠️  Warning: There are uncommitted changes")
        print("Uncommitted files:")
        for line in status.split('\n'):
            if line.strip():
                print(f"   {line}")
        return False
    
    return True

def sync_all_versions():
    """Sync all version references"""
    print(f"🔄 Syncing all version references to {__version__}...")
    
    # Update all versions using version manager
    result = run_command(f"{PYTHON} scripts/version_manager.py set {__version__}")
    if result is None:
        print("❌ Failed to update versions")
        return False
    
    print("✅ All version references synchronized")
    return True

def check_version_consistency():
    """Check if all version references are consistent"""
    print("🔍 Checking version consistency...")
    
    issues = []
    
    # Check README.md
    readme_file = PROJECT_ROOT / "README.md"
    if readme_file.exists():
        with open(readme_file, 'r') as f:
            content = f.read()
            if f"v{__version__}" not in content and f"Version {__version__}" not in content:
                issues.append("README.md version reference")
    
    # Check predict.py
    predict_file = PROJECT_ROOT / "rxnrecer/cli/predict.py"
    if predict_file.exists():
        with open(predict_file, 'r') as f:
            content = f.read()
            if f"v{__version__}" not in content:
                issues.append("predict.py version reference")
    
    # Check docs files
    docs_files = ["docs/INSTALL.md", "docs/RELEASE_NOTES.md"]
    for doc_file in docs_files:
        file_path = PROJECT_ROOT / doc_file
        if file_path.exists():
            with open(file_path, 'r') as f:
                content = f.read()
                if f"v{__version__}" not in content and f"Version {__version__}" not in content:
                    issues.append(f"{doc_file} version reference")
    
    if issues:
        print("⚠️  Version consistency issues found:")
        for issue in issues:
            print(f"   - {issue}")
        return False
    else:
        print("✅ All version references are consistent")
        return True

def main():
    """Main function"""
    print("🔄 RXNRECer Version Synchronization")
    print(f"Current version: {__version__}")
    print("=" * 50)
    
    # Check git status
    if not check_git_status():
        print("💡 Please commit or stash your changes first")
        sys.exit(1)
    
    # Sync all versions
    if not sync_all_versions():
        print("❌ Version synchronization failed")
        sys.exit(1)
    
    # Check consistency
    if not check_version_consistency():
        print("❌ Version consistency check failed")
        sys.exit(1)
    
    print("\n✅ Version synchronization completed successfully!")
    print("💡 You can now commit your changes:")
    print("   git add .")
    print("   git commit -m 'chore: Sync version references'")

if __name__ == '__main__':
    main()

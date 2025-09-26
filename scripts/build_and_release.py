#!/usr/bin/env python3
"""
Build and release script for RXNRECer package.
This script automates the process of building and releasing the package.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

def run_command(command, check=True, capture_output=False):
    """Run a shell command."""
    print(f"Running: {command}")
    try:
        result = subprocess.run(
            command, 
            shell=True, 
            check=check, 
            capture_output=capture_output,
            text=True
        )
        if capture_output:
            return result.stdout.strip()
        return result
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")
        if check:
            sys.exit(1)
        return None

def clean_build():
    """Clean previous build artifacts."""
    print("🧹 Cleaning previous build artifacts...")
    
    dirs_to_clean = ["build", "dist", "*.egg-info"]
    for pattern in dirs_to_clean:
        if "*" in pattern:
            # Handle wildcard patterns
            for item in Path(".").glob(pattern):
                if item.is_dir():
                    shutil.rmtree(item)
                else:
                    item.unlink()
        else:
            if Path(pattern).exists():
                if Path(pattern).is_dir():
                    shutil.rmtree(pattern)
                else:
                    Path(pattern).unlink()
    
    print("✅ Build artifacts cleaned")

def install_build_deps():
    """Install build dependencies."""
    print("📦 Installing build dependencies...")
    
    build_deps = [
        "build",
        "wheel", 
        "setuptools",
        "setuptools_scm[toml]",
        "twine"
    ]
    
    for dep in build_deps:
        run_command(f"pip install {dep}")
    
    print("✅ Build dependencies installed")

def build_package():
    """Build the package."""
    print("🔨 Building package...")
    
    # Build source distribution and wheel
    run_command("python -m build")
    
    # Check if build was successful
    if not Path("dist").exists():
        print("❌ Build failed - dist directory not created")
        sys.exit(1)
    
    print("✅ Package built successfully")

def check_package():
    """Check the built package."""
    print("🔍 Checking package...")
    
    # Check source distribution
    sdist_files = list(Path("dist").glob("*.tar.gz"))
    if sdist_files:
        run_command(f"twine check {sdist_files[0]}")
    
    # Check wheel
    wheel_files = list(Path("dist").glob("*.whl"))
    if wheel_files:
        run_command(f"twine check {wheel_files[0]}")
    
    print("✅ Package check completed")

def test_installation():
    """Test package installation."""
    print("🧪 Testing package installation...")
    
    # Create test environment
    test_dir = Path("test_install")
    if test_dir.exists():
        shutil.rmtree(test_dir)
    test_dir.mkdir()
    
    # Copy wheel to test directory
    wheel_files = list(Path("dist").glob("*.whl"))
    if wheel_files:
        shutil.copy(wheel_files[0], test_dir)
        
        # Test install
        os.chdir(test_dir)
        run_command(f"pip install {wheel_files[0].name}")
        
        # Test import
        run_command("python -c 'import rxnrecer; print(f\"RXNRECer {rxnrecer.__version__} installed successfully\")'")
        
        # Test CLI
        run_command("rxnrecer --help")
        
        os.chdir("..")
        shutil.rmtree(test_dir)
        print("✅ Package installation test passed")
    else:
        print("❌ No wheel file found for testing")

def upload_to_test_pypi():
    """Upload package to Test PyPI."""
    print("🚀 Uploading to Test PyPI...")
    
    # Check if TWINE_USERNAME and TWINE_PASSWORD are set
    if not (os.getenv("TWINE_USERNAME") and os.getenv("TWINE_PASSWORD")):
        print("⚠️  TWINE_USERNAME and TWINE_PASSWORD not set")
        print("   Set these environment variables to upload to Test PyPI")
        return False
    
    try:
        run_command("twine upload --repository testpypi dist/*")
        print("✅ Package uploaded to Test PyPI successfully")
        return True
    except Exception as e:
        print(f"❌ Upload to Test PyPI failed: {e}")
        return False

def upload_to_pypi():
    """Upload package to PyPI."""
    print("🚀 Uploading to PyPI...")
    
    # Check if TWINE_USERNAME and TWINE_PASSWORD are set
    if not (os.getenv("TWINE_USERNAME") and os.getenv("TWINE_PASSWORD")):
        print("⚠️  TWINE_USERNAME and TWINE_PASSWORD not set")
        print("   Set these environment variables to upload to PyPI")
        return False
    
    # Confirm upload
    response = input("Are you sure you want to upload to PyPI? (y/N): ")
    if response.lower() != 'y':
        print("Upload cancelled")
        return False
    
    try:
        run_command("twine upload dist/*")
        print("✅ Package uploaded to PyPI successfully")
        return True
    except Exception as e:
        print(f"❌ Upload to PyPI failed: {e}")
        return False

def create_release_notes():
    """Create release notes."""
    print("📝 Creating release notes...")
    
    # Get version from package
    try:
        version = run_command("python -c 'import rxnrecer; print(rxnrecer.__version__)'", capture_output=True)
        print(f"Package version: {version}")
    except:
        version = "1.3.0"
    
    release_notes = f"""# RXNRECer {version} Release Notes

## 🎉 New Release

This release includes:
- Complete package structure for pip installation
- Command-line interface (CLI) for easy use
- Python API for programmatic access
- Comprehensive documentation and examples

## 📦 Installation

```bash
# Install from PyPI
pip install rxnrecer

# Install with development dependencies
pip install rxnrecer[dev]

# Install from GitHub
pip install git+https://github.com/kingstdio/RXNRECer.git
```

## 🚀 Usage

```bash
# Basic prediction
rxnrecer -i input.fasta -o output.tsv -m s1

# Get help
rxnrecer --help
```

## 📚 Documentation

- [Installation Guide](docs/INSTALL.md)
- [Release Notes](docs/RELEASE_NOTES.md)

## 🔗 Links

- **PyPI**: https://pypi.org/project/rxnrecer/
- **GitHub**: https://github.com/kingstdio/RXNRECer
- **Documentation**: https://github.com/kingstdio/RXNRECer#readme

---
Released on: {os.popen('date').read().strip()}
"""
    
    with open("RELEASE_NOTES_PYPI.md", "w") as f:
        f.write(release_notes)
    
    print("✅ Release notes created: RELEASE_NOTES_PYPI.md")

def main():
    """Main function."""
    print("🚀 RXNRECer Build and Release Script")
    print("=" * 50)
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        command = sys.argv[1]
    else:
        command = "build"
    
    if command == "clean":
        clean_build()
    elif command == "build":
        clean_build()
        install_build_deps()
        build_package()
        check_package()
        test_installation()
    elif command == "test":
        test_installation()
    elif command == "upload-test":
        build_package()
        upload_to_test_pypi()
    elif command == "upload":
        build_package()
        upload_to_pypi()
    elif command == "release":
        clean_build()
        install_build_deps()
        build_package()
        check_package()
        test_installation()
        create_release_notes()
        print("\n🎉 Release preparation completed!")
        print("Next steps:")
        print("1. Review the built package in dist/")
        print("2. Test installation: python scripts/build_and_release.py test")
        print("3. Upload to Test PyPI: python scripts/build_and_release.py upload-test")
        print("4. Upload to PyPI: python scripts/build_and_release.py upload")
    else:
        print(f"Unknown command: {command}")
        print("Available commands:")
        print("  clean      - Clean build artifacts")
        print("  build      - Build package")
        print("  test       - Test package installation")
        print("  upload-test - Upload to Test PyPI")
        print("  upload     - Upload to PyPI")
        print("  release    - Full release preparation")

if __name__ == "__main__":
    main()

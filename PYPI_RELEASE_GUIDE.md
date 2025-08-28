# RXNRECer PyPI Release Guide

This guide explains how to release RXNRECer to PyPI for easy pip installation.

## 🎯 Overview

RXNRECer is now configured as a proper Python package that can be:
- **Installed via pip**: `pip install rxnrecer`
- **Installed from GitHub**: `pip install git+https://github.com/kingstdio/RXNRECer.git`
- **Built and distributed**: Source and wheel distributions

## 📦 Package Structure

### Core Package Files
- `setup.py` - Traditional setup configuration
- `pyproject.toml` - Modern Python packaging configuration
- `requirements.txt` - Dependencies with version constraints
- `MANIFEST.in` - Files to include in distribution

### Package Configuration
- **Name**: `rxnrecer`
- **Version**: `1.0.0`
- **Author**: Zhenkun Shi
- **License**: MIT
- **Python**: 3.10+
- **Entry Point**: `rxnrecer` command

## 🚀 Quick Release

### 1. Build Package
```bash
# Build source and wheel distributions
python scripts/build_and_release.py build

# This will:
# - Clean previous builds
# - Install build dependencies
# - Build source distribution (.tar.gz)
# - Build wheel distribution (.whl)
# - Check package integrity
# - Test installation
```

### 2. Test Package
```bash
# Test the built package
python scripts/build_and_release.py test

# This will:
# - Install the wheel in a test environment
# - Verify import works
# - Test CLI command
# - Clean up test environment
```

### 3. Prepare for Release
```bash
# Complete release preparation
python scripts/build_and_release.py release

# This will:
# - Build package
# - Run all tests
# - Create release notes
# - Provide next steps
```

## 📤 Publishing to PyPI

### Prerequisites
1. **PyPI Account**: Create account at https://pypi.org/
2. **Test PyPI Account**: Create account at https://test.pypi.org/
3. **API Tokens**: Generate API tokens for both accounts

### Environment Setup
```bash
# Set PyPI credentials
export TWINE_USERNAME="your_pypi_username"
export TWINE_PASSWORD="your_pypi_api_token"

# Or create .pypirc file
cat > ~/.pypirc << EOF
[pypi]
username = your_pypi_username
password = your_pypi_api_token

[testpypi]
username = your_pypi_username
password = your_test_pypi_api_token
EOF
```

### Upload to Test PyPI
```bash
# Upload to Test PyPI first
python scripts/build_and_release.py upload-test

# This will:
# - Build package
# - Upload to Test PyPI
# - Verify upload success
```

### Upload to PyPI
```bash
# Upload to production PyPI
python scripts/build_and_release.py upload

# This will:
# - Build package
# - Prompt for confirmation
# - Upload to PyPI
# - Verify upload success
```

## 🔧 Manual Build and Upload

### Build Manually
```bash
# Install build tools
pip install build wheel twine

# Build package
python -m build

# Check package
twine check dist/*
```

### Upload Manually
```bash
# Upload to Test PyPI
twine upload --repository testpypi dist/*

# Upload to PyPI
twine upload dist/*
```

## 📋 Release Checklist

### Before Release
- [ ] Update version in `rxnrecer/__init__.py`
- [ ] Update version in `setup.py`
- [ ] Update version in `pyproject.toml`
- [ ] Test package builds successfully
- [ ] Test package installs correctly
- [ ] Test CLI works as expected
- [ ] Verify all dependencies are correct

### During Release
- [ ] Build package: `python scripts/build_and_release.py build`
- [ ] Test package: `python scripts/build_and_release.py test`
- [ ] Upload to Test PyPI: `python scripts/build_and_release.py upload-test`
- [ ] Test installation from Test PyPI
- [ ] Upload to PyPI: `python scripts/build_and_release.py upload`

### After Release
- [ ] Verify package appears on PyPI
- [ ] Test installation: `pip install rxnrecer`
- [ ] Test CLI: `rxnrecer --help`
- [ ] Update GitHub release notes
- [ ] Announce release to users

## 🐛 Troubleshooting

### Build Issues
```bash
# Clean and rebuild
python scripts/build_and_release.py clean
python scripts/build_and_release.py build
```

### Upload Issues
```bash
# Check credentials
echo $TWINE_USERNAME
echo $TWINE_PASSWORD

# Test connection
twine check dist/*
```

### Installation Issues
```bash
# Test in clean environment
python -m venv test_env
source test_env/bin/activate
pip install rxnrecer
rxnrecer --help
```

## 📊 Package Information

### Distribution Files
- **Source Distribution**: `rxnrecer-1.0.0.tar.gz` (~79KB)
- **Wheel Distribution**: `rxnrecer-1.0.0-py3-none-any.whl` (~65KB)

### Dependencies
- **Core**: PyTorch, NumPy, Pandas, Scikit-learn
- **ML**: Transformers, Fair-ESM
- **Bio**: Biopython, RDKit
- **Utils**: Click, Tqdm, PyYAML

### Optional Dependencies
- **Dev**: pytest, black, flake8, isort
- **Docs**: sphinx, sphinx-rtd-theme
- **Full**: All optional dependencies

## 🔗 Useful Links

- **PyPI**: https://pypi.org/project/rxnrecer/
- **Test PyPI**: https://test.pypi.org/project/rxnrecer/
- **GitHub**: https://github.com/kingstdio/RXNRECer
- **Documentation**: https://github.com/kingstdio/RXNRECer#readme

## 📞 Support

If you encounter issues during release:
1. Check the build logs for errors
2. Verify all dependencies are available
3. Test in a clean environment
4. Contact: zhenkun.shi@tib.cas.cn

---

**Happy releasing! 🚀**

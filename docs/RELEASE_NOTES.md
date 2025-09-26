# RXNRECer Release Notes

## Version 1.3.5 - PyPI Release

**Release Date**: August 2025  
**Status**: Production Ready

### 🎉 Major Features

- **PyPI Package**: Now available as `pip install rxnrecer`
- **Smart Caching**: Automatic result caching for faster repeated predictions
- **Multi-Stage Prediction**: S1 (reaction), S2 (integration), S3 (LLM reasoning) modes
- **ESM-2 Integration**: State-of-the-art protein language model embeddings
- **GPU Acceleration**: Full CUDA support for high-performance inference

### 🚀 New Capabilities

#### PyPI Installation
```bash
pip install rxnrecer
```
- ✅ One-command installation
- ✅ Automatic dependency management
- ✅ Global availability
- ✅ Easy updates

#### Smart Caching System
- **Automatic caching**: Results cached by input file, mode, and format
- **Cache management**: `rxnrecer-cache` command for status and cleanup
- **Performance boost**: Instant results for repeated predictions

#### Command Line Interface
- `rxnrecer` - Main prediction command
- `rxnrecer-download-data` - Data and model download
- `rxnrecer-cache` - Cache management

### 📁 Project Structure

```
RXNRECer/                    # Project root
├── rxnrecer/               # Main Python package
│   ├── cli/                # Command-line interface
│   ├── config/             # Configuration
│   ├── lib/                # Core libraries
│   └── utils/              # Utility functions
├── data/                    # Data files (download required)
├── ckpt/                   # Model checkpoints (download required)
├── results/                 # Output results
├── docs/                    # Documentation
└── scripts/                 # Build scripts
```

### 🔧 Technical Requirements

- Python 3.10+
- PyTorch 2.0+
- CUDA 11.0+ (recommended)
- 32GB+ RAM
- 40GB+ disk space

### 🆕 New Features

#### 📥 Automatic Data Download
```bash
rxnrecer-download-data        # Download all files (~35.8GB)
rxnrecer-download-data --data-only      # Data only (~8.8GB)
rxnrecer-download-data --models-only    # Models only (~14GB)
rxnrecer-download-data --extools-only   # External tools only (~13GB)
```

#### 💾 Smart Caching
```bash
rxnrecer-cache status         # Check cache status
rxnrecer-cache clear --all    # Clear all cache
```

### 🐛 Bug Fixes

- Fixed LLM API configuration variables
- Improved error handling in cache system
- Cleaned up unused configuration variables
- Fixed path resolution issues

### 📚 Documentation

- Updated README for PyPI installation
- Simplified installation guide
- Added cache management examples
- Streamlined project structure

### 🔄 Migration from v1.3.5

1. **Installation**: Use `pip install rxnrecer` instead of GitHub clone
2. **Caching**: New cache system automatically improves performance
3. **Configuration**: Simplified LLM API setup
4. **CLI**: Enhanced commands with better help messages

### 🚀 Quick Start

```bash
# Install
pip install rxnrecer

# Download data
rxnrecer-download-data

# Run prediction
rxnrecer -i input.fasta -o output.tsv -m s1
```

### 📞 Support

- **PyPI**: https://pypi.org/project/rxnrecer/
- **Documentation**: https://github.com/kingstdio/RXNRECer#readme
- **Contact**: zhenkun.shi@tib.cas.cn

---

**Get started: `pip install rxnrecer`**

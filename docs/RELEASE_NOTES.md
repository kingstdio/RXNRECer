# RXNRECer Release Notes

## Version 1.4.2 - Public Package Cleanup

**Release Date**: July 21, 2026  
**Status**: Production Ready

### Updates

- Cleaned the source distribution to include only user-facing documentation.
- Removed non-user-facing workflow notes and development templates from the public release tree.
- Simplified release notes to focus on installation and runtime changes relevant to end users.

---

## Version 1.4.1 - Installation Compatibility

**Release Date**: July 20, 2026  
**Status**: Production Ready

### Fixes

- Constrained package metadata to `numpy>=1.21.0,<2.0` to match the validated environment.
- Improved model checkpoint compatibility for released runtime assets.

---

## Version 1.4.0 - Documentation and Packaging Refresh

**Release Date**: July 20, 2026  
**Status**: Production Ready

### Major Updates

- Added clearer installation and reproducibility guidance.
- Documented the validated environment profile, including Python, NumPy, PyTorch, CUDA, pandas, scikit-learn, and Apptainer versions.
- Clarified S2/S3 container-runtime requirements for CatFam and ECRECer `.sif` images.
- Added public webserver documentation and local usage guidance.
- Added LLM API setup documentation for hosted providers and local vLLM/SGLang-compatible services.
- Restored dynamic package version generation from `rxnrecer.__version__` for release builds.

### Migration from v1.3.7

1. Install or upgrade with `python -m pip install -U rxnrecer`.
2. Use NumPy 1.x for the current validated release profile.
3. For S2/S3, verify `singularity --version` and inspect the downloaded CatFam/ECRECer images before large runs.

---

## Version 1.3.7 - SVG Generation & Stability Improvements

**Release Date**: November 2025  
**Status**: Production Ready

### 🎉 Major Features

- **New CLI Tool**: `rxnrecer-genmolsvg` - Batch generate molecule SVGs for all reactions
- **Deterministic SVG Filenames**: Stable file naming based on CHEBI ID or canonical SMILES
- **SVG Path Fix**: Fixed missing SVG files in JSON output
- **Improved Reliability**: Cache disabled by default to ensure fresh SVG generation

### 🚀 New Capabilities

#### Molecule SVG Generation
```bash
rxnrecer-genmolsvg              # Generate all molecule SVGs
rxnrecer-genmolsvg --limit 100  # Test with first 100 reactions
```

#### Stable File Naming
- SVG files named using CHEBI ID (if available) or canonical SMILES hash
- Same compound always generates same filename
- Compatible with pre-generated SVG batches

### 🐛 Bug Fixes

- Fixed SVG file path resolution in JSON output
- Fixed missing SVG files in reaction details
- Improved project root directory detection
- Better error handling for molecule generation

### 📚 Documentation

- Updated version references to 1.3.7
- Added SVG generation documentation
- Improved CLI help messages

### 🔄 Migration from v1.3.6

1. **SVG Files**: Run `rxnrecer-genmolsvg` to pre-generate all molecule SVGs
2. **Cache**: Cache is now disabled by default for reliability
3. **CLI**: New `rxnrecer-genmolsvg` command available

---

## Version 1.3.7 - PyPI Release

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
└── docs/                    # User documentation
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

### 🔄 Migration from v1.3.6

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

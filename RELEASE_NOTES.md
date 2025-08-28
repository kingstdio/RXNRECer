# RXNRECer Release Notes

## Version 1.0.0 - Production Release

**Release Date**: January 2025  
**Branch**: `release` (Default branch)

### 🎉 Major Features

- **Production-Ready CLI**: Complete command-line interface for easy deployment
- **Multi-Stage Prediction**: S1 (reaction), S2 (integration), S3 (LLM reasoning) modes
- **ESM-2 Integration**: State-of-the-art protein language model embeddings
- **GPU Acceleration**: Full CUDA support for high-performance inference
- **Comprehensive Documentation**: Complete user guide and API documentation

### 🚀 New Capabilities

#### Command Line Interface
- `rxnrecer` command with intuitive options
- Support for FASTA input files
- Multiple output formats (TSV, JSON)
- Batch processing with configurable batch sizes

#### Prediction Modes
- **S1 Mode**: Basic reaction prediction using ESM-2
- **S2 Mode**: Enhanced prediction with reaction details
- **S3 Mode**: Advanced prediction with LLM reasoning

#### Model Architecture
- ESM-2 protein sequence embeddings
- Advanced neural network architectures
- Optimized for enzyme reaction prediction
- Support for large protein datasets

### 📁 Project Structure

```
rxnrecer/
├── cli/                    # Command-line interface
├── config/                 # Configuration management
├── lib/                    # Core library modules
├── models/                 # Neural network models
├── preprocessing/           # Data preprocessing
├── utils/                  # Utility functions
├── data/                   # Data directory
└── ckpt/                   # Model checkpoints
```

### 🔧 Technical Requirements

- Python 3.10+
- PyTorch 2.0+
- CUDA 11.0+ (recommended)
- 32GB+ RAM
- 40GB+ disk space

### 📊 Performance Improvements

- Optimized memory usage for large datasets
- Efficient batch processing
- GPU acceleration support
- Streamlined data pipeline

### 🐛 Bug Fixes

- Removed deprecated training scripts
- Cleaned up experimental notebooks
- Fixed configuration issues
- Improved error handling

### 📚 Documentation

- Comprehensive README with examples
- Installation and setup guide
- API documentation
- Usage examples for all modes

### 🔄 Migration Guide

This is a major release with significant changes:

1. **New CLI**: Use `rxnrecer` command instead of Python scripts
2. **Simplified Structure**: Removed experimental and training code
3. **Production Focus**: Optimized for deployment and inference
4. **Better Configuration**: Centralized configuration management

### 📥 Installation

```bash
# Clone the repository
git clone https://github.com/kingstdio/RXNRECer.git
cd RXNRECer

# Checkout release branch
git checkout release

# Install in development mode
pip install -e .
```

### 🚀 Quick Start

```bash
# Basic prediction
rxnrecer -i input.fasta -o output.tsv -m s1

# Advanced prediction with LLM
rxnrecer -i input.fasta -o output.json -m s3 -f json
```

### 🤝 Contributing

- Fork the repository
- Create feature branches from `release`
- Submit pull requests
- Follow coding standards

### 📞 Support

- **Issues**: [GitHub Issues](https://github.com/kingstdio/RXNRECer/issues)
- **Documentation**: [README.md](README.md)
- **Contact**: zhenkun.shi@tib.cas.cn

---

**Note**: This release branch is now the default branch and contains the production-ready version of RXNRECer.

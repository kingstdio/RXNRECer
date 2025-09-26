# RXNRECer

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/badge/PyPI-rxnrecer-blue.svg)](https://pypi.org/project/rxnrecer/)

**RXNRECer {{VERSION}}** is a deep learning framework for predicting enzyme-catalyzed reactions from protein sequences.
It is the official implementation of "RXNRECer: Active Learning with Protein Language Models for Fine-Grained Enzyme Reaction Prediction."

**🎉 Now available on PyPI for easy installation!**

## 🚀 Features

- **Multi-Stage Prediction**: S1 (reaction prediction), S2 (reaction integration), S3 (LLM reasoning)
- **Protein Sequence Analysis**: Process protein sequences in FASTA format
- **Deep Learning Models**: ESM-2 embeddings with advanced neural architectures
- **GPU Acceleration**: CUDA support for faster inference
- **Easy-to-use CLI**: Simple command-line interface with comprehensive options
- **Smart Caching**: Automatic result caching for faster repeated predictions

## 📋 Requirements

- Python 3.10+
- PyTorch 2.0+
- CUDA 11.0+ (recommended)
- 32GB+ RAM
- 40GB+ disk space

## 🚀 Quick Start

### 1. Install (Recommended)

```bash
# Install from PyPI (recommended)
pip install rxnrecer=={{VERSION}}

# Or install from GitHub
pip install git+https://github.com/kingstdio/RXNRECer.git
```

### 2. Download Data

```bash
# Download required data and model files (~35.8GB total)
rxnrecer-download-data

# Or download separately
rxnrecer-download-data --data-only      # ~8.8GB
rxnrecer-download-data --models-only    # ~14GB
rxnrecer-download-data --extools-only   # ~13GB
```

### 3. Run Prediction

```bash
# Basic prediction (S1 mode)
rxnrecer -i input.fasta -o output.tsv -m s1

# Detailed prediction (S2 mode)
rxnrecer -i input.fasta -o output.tsv -m s2

# LLM reasoning (S3 mode, requires API key)
rxnrecer -i input.fasta -o output.json -m s3 -f json
```

## 📦 Installation Options

### PyPI Installation (Recommended)
```bash
pip install rxnrecer=={{VERSION}}
```

### GitHub Installation (Latest)
```bash
pip install git+https://github.com/kingstdio/RXNRECer.git
```

## 📚 Documentation

- **[Installation Guide](docs/INSTALL.md)** - Detailed setup instructions
- **[Release Notes](docs/RELEASE_NOTES.md)** - Version information

## 📞 Contact

- **Author**: Zhenkun Shi
- **Email**: zhenkun.shi@tib.cas.cn
- **Project**: [https://github.com/kingstdio/RXNRECer](https://github.com/kingstdio/RXNRECer)
- **PyPI**: [https://pypi.org/project/rxnrecer/](https://pypi.org/project/rxnrecer/)

---

**🎯 Get started now with: `pip install rxnrecer=={{VERSION}}`**

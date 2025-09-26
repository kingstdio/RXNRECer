# RXNRECer

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/badge/PyPI-rxnrecer-blue.svg)](https://pypi.org/project/rxnrecer/)

**RXNRECer v1.3.4** is a deep learning framework for predicting enzyme-catalyzed reactions from protein sequences.
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
pip install rxnrecer

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

## 🔧 Usage

### Command Line Options

```bash
rxnrecer [OPTIONS]

Options:
  -i, --input_fasta    Input FASTA file path (required)
  -o, --output_file    Output file path
  -f, --format         Output format: tsv or json (default: tsv)
  -m, --mode           Prediction mode: s1, s2, or s3 (default: s1)
  -b, --batch_size     Batch size for processing (default: 100)
  -v, --version        Show version
```

### Examples

```bash
# Basic usage
rxnrecer -i proteins.fasta -o results.tsv

# Custom batch size
rxnrecer -i proteins.fasta -o results.tsv -b 50

# JSON output
rxnrecer -i proteins.fasta -o results.json -f json

# Use default output path
rxnrecer -i proteins.fasta -m s1
```

### Input Format

FASTA file with protein sequences:

```
>P12345|Sample protein 1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>P67890|Sample protein 2
MKLIVWALLLLAAWAVERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
```

### Output Formats

**TSV Output (S1/S2):**
```tsv
input_id	RXNRECer	RXNRECer_with_prob	rxn_details
P12345	RHEA:24076;RHEA:14709	0.9999;0.9999	[reaction details]
```

**JSON Output (S3):**
```json
[
  {
    "reaction_id": "RHEA:24076",
    "prediction_confidence": 0.9999,
    "reaction_details": {...}
  }
]
```

## 🆕 Advanced Features

### Smart Caching
Results are automatically cached for faster repeated predictions:

```bash
# Check cache status
rxnrecer-cache status

# Clear cache
rxnrecer-cache clear --all
```

### Data Management
Easy data and model file management:

```bash
# Download data
rxnrecer-download-data

# Force re-download
rxnrecer-download-data --force
```

## 📁 Project Structure

```
RXNRECer/                               # Project root (release)
├── .github/                            # CI/CD workflows
│   └── workflows/
├── rxnrecer/                           # Main Python package
│   ├── cli/                            # Command-line interface
│   ├── config/                         # Configuration
│   ├── lib/                            # Core libraries
│   │   ├── datasource/                 # Data source handling
│   │   ├── embedding/                  # Protein embeddings
│   │   ├── evaluation/                 # Evaluation helpers
│   │   ├── llm/                        # Language model integration
│   │   ├── ml/                         # Machine learning utilities
│   │   ├── model/                      # Model architectures
│   │   ├── rxn/                        # Reaction processing
│   │   └── smi/                        # SMILES handling
│   ├── models/                         # Model wrappers
│   └── utils/                          # Utility functions
│
├── extools/                            # External tools (downloaded)
│   ├── ec/                             # EC-related resources
│   └── msa/                            # MSA binaries (e.g., diamond)
│
├── data/                               # Data files (download required)
│   ├── chebi/                          # ChEBI database
│   ├── cpd_svg/                        # Compound SVG files
│   ├── datasets/                       # Training datasets
│   ├── dict/                           # Dictionary files
│   ├── feature_bank/                   # Feature bank
│   ├── rhea/                           # RHEA database
│   ├── rxn_json/                       # Reaction JSON files
│   ├── sample/                         # Sample data
│   └── uniprot/                        # UniProt database
│
├── ckpt/                              # Model checkpoints (download required)
│   ├── esm/                           # ESM models
│   ├── prostt5/                       # ProSTT5 models
│   └── rxnrecer/                      # RXNRECer model files
│
├── results/                            # Output results
│   ├── cache/                          # Prediction cache
│   ├── logs/                           # Log files
│   ├── predictions/                    # Prediction outputs
│   └── sample/                         # Sample results
│
├── docs/                               # Documentation
├── scripts/                            # Build and utility scripts
├── MANIFEST.in                         # Package data manifest
├── pyproject.toml                      # Build and dependencies for PyPI
├── environment_rxnrecer-release.yml    # Conda environment
├── LICENSE                             # MIT License
├── README.md                           # This file
└── .gitignore                          # Git ignore rules
```

## 🔧 Configuration

For S3 mode (LLM reasoning), set your API key:

```bash
export LLM_API_KEY="your_api_key_here"
export LLM_API_URL="your_api_url_here"
```

**Examples:**

```bash
# OpenRouter
export LLM_API_KEY="sk-or-v1-your_openrouter_key_here"
export LLM_API_URL="https://openrouter.ai/api/v1"

# OpenAI
export LLM_API_KEY="sk-your_openai_key_here"
export LLM_API_URL="https://api.openai.com/v1"

# Anthropic
export LLM_API_KEY="sk-ant-your_anthropic_key_here"
export LLM_API_URL="https://api.anthropic.com"
```

### Jupyter Notebook Setup

```python
import os
from rxnrecer.config import config as cfg

# Set your API credentials
cfg.LLM_API_KEY = "your_api_key_here"
cfg.LLM_API_URL = "your_api_url_here"
```

## 📦 Installation Options

### PyPI Installation (Recommended)
```bash
pip install rxnrecer
```

### GitHub Installation (Latest)
```bash
pip install git+https://github.com/kingstdio/RXNRECer.git
```
- 🔧 **Development**: Latest development version
- 🔧 **Custom**: For advanced users

## 📚 Documentation

- **[Installation Guide](docs/INSTALL.md)** - Detailed setup instructions
- **[Release Notes](docs/RELEASE_NOTES.md)** - Version information

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Open a Pull Request

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 📞 Contact

- **Author**: Zhenkun Shi
- **Email**: zhenkun.shi@tib.cas.cn
- **Project**: [https://github.com/kingstdio/RXNRECer](https://github.com/kingstdio/RXNRECer)
- **PyPI**: [https://pypi.org/project/rxnrecer/](https://pypi.org/project/rxnrecer/)

---

**🎯 Get started now with: `pip install rxnrecer`**





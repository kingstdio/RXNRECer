# RXNRECer

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**RXNRECer** is a deep learning framework for predicting enzyme-catalyzed reactions from protein sequences.
It is the official implementation of "RXNRECer: Active Learning with Protein Language Models for Fine-Grained Enzyme Reaction Prediction."

## 🚀 Features

- **Multi-Stage Prediction**: S1 (reaction prediction), S2 (reaction integration), S3 (LLM reasoning)
- **Protein Sequence Analysis**: Process protein sequences in FASTA format
- **Deep Learning Models**: ESM-2 embeddings with advanced neural architectures
- **GPU Acceleration**: CUDA support for faster inference
- **Easy-to-use CLI**: Simple command-line interface with comprehensive options

## 📋 Requirements

- Python 3.10+
- PyTorch 2.0+
- CUDA 11.0+ (recommended for GPU acceleration)
- 32GB+ RAM (for large protein datasets)
- 40GB+ disk space (for models and data)

## 🚀 Quick Start

### 1. Install RXNRECer

#### Option 1: Install from PyPI (Recommended)

```bash
# Install the latest stable version
pip install rxnrecer

# Install with development dependencies
pip install rxnrecer[dev]

# Install with all dependencies
pip install rxnrecer[full]
```

#### Option 2: Install from GitHub

```bash
# Install directly from GitHub
pip install git+https://github.com/kingstdio/RXNRECer.git

# Install specific branch
pip install git+https://github.com/kingstdio/RXNRECer.git@release

# Install specific version
pip install git+https://github.com/kingstdio/RXNRECer.git@v1.0.0
```

#### Option 3: Install in Development Mode

```bash
# Clone the repository
git clone https://github.com/kingstdio/rxnrecer.git
cd rxnrecer

# Install in development mode
pip install -e .
```

### 2. Download Data and Model Files

**Important**: Due to GitHub file size limits, you need to download data and model files separately.

```bash
# Download data files (~8.6GB)
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz"

# Download model files (~11.9GB)  
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz"

# Extract to current directory (data.tar.gz -> data/, ckpt.tar.gz -> ckpt/)
tar -xzf data.tar.gz
tar -xzf ckpt.tar.gz
```

### 3. Configure LLM API (Required for S3 mode)

```bash
# Set your API key
export LLM_API_KEY="your_api_key_here"

# Optional: Set custom API endpoint
export LLM_API_URL="https://openrouter.ai/api/v1"
```

If you're using Jupyter notebooks, you need to set environment variables within the notebook:

```python
# At the beginning of your notebook, run this:
from rxnrecer.config import config as cfg
cfg.LLM_API_KEY = "your_actual_api_key_here"
cfg.LLM_API_URL = "https://openrouter.ai/api/v1"
```

### 4. Run Prediction

```bash
# Basic S1 prediction
rxnrecer -i data/sample/sample10.fasta -o results/output.tsv -m s1

# S2 prediction with reaction details
rxnrecer -i data/sample/sample10.fasta -o results/output.tsv -m s2

# S3 prediction with LLM reasoning
rxnrecer -i data/sample/sample10.fasta -o results/output.json -m s3 -f json
```

## 📁 Project Structure

```
rxnrecer/                          # Main Python package
├── cli/                           # Command-line interface
├── config/                        # Configuration files
├── lib/                           # Core library modules
├── models/                        # Neural network models
├── preprocessing/                  # Data preprocessing
├── utils/                         # Utility functions
├── data/                          # Data directory (download from S3)
└── ckpt/                          # Model checkpoints (download from S3)
```

## 🔧 Configuration

### Model Modes

- **S1 Mode**: Basic reaction prediction using ESM-2 embeddings
- **S2 Mode**: Enhanced prediction with reaction details and equations
- **S3 Mode**: Advanced prediction with LLM reasoning and confidence scoring

### Command Line Options

```bash
rxnrecer [OPTIONS]

Options:
  -i, --input_fasta    Input FASTA file path
  -o, --output_file    Output file path
  -f, --format         Output format: tsv or json
  -m, --mode           Prediction mode: s1, s2, or s3
  -b, --batch_size     Batch size for processing (default: 100)
```

## 📊 Input/Output Format

### FASTA Input

```
>P12345|Sample protein 1
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
```

### TSV Output (S1/S2 Mode)

```tsv
input_id	RXNRECer	RXNRECer_with_prob	rxn_details
P12345	RHEA:24076;RHEA:14709	0.9999;0.9999	[reaction details]
```

### JSON Output (S3 Mode)

```json
[
  {
    "reaction_id": "RHEA:24076",
    "prediction_confidence": 0.9999,
    "reaction_details": {...},
    "reaction_rxnrecer_s3": {...}
  }
]
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📦 Installation Options

### PyPI Installation (Recommended)
```bash
pip install rxnrecer
```

### GitHub Installation
```bash
pip install git+https://github.com/kingstdio/RXNRECer.git
```

### Development Installation
```bash
git clone https://github.com/kingstdio/RXNRECer.git
cd RXNRECer
pip install -e .
```

## 📚 Additional Documentation

- **[Installation Guide](INSTALL.md)** - Detailed installation instructions
- **[Release Notes](RELEASE_NOTES.md)** - Version information and changes
- **[Deployment Guide](DEPLOYMENT_GUIDE.md)** - Production deployment guide

## 🔧 Build and Release

To build and release the package:

```bash
# Build package
python scripts/build_and_release.py build

# Test installation
python scripts/build_and_release.py test

# Prepare for release
python scripts/build_and_release.py release
```

## 📞 Contact

- **Author**: Zhenkun Shi
- **Email**: zhenkun.shi@tib.cas.cn
- **Project**: [https://github.com/kingstdio/rxnrecer](https://github.com/kingstdio/rxnrecer)
- **PyPI**: [https://pypi.org/project/rxnrecer/](https://pypi.org/project/rxnrecer/)





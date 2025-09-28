# RXNRECer Installation Guide

**Version 1.3.6**

## 🚀 Quick Installation

### Install from PyPI (Recommended)

```bash
pip install rxnrecer
```

### Install from GitHub

```bash
pip install git+https://github.com/kingstdio/RXNRECer.git
```

## 📋 Prerequisites

- **Python**: 3.10+
- **RAM**: 32GB+ (recommended)
- **Disk Space**: 40GB+ (for models and data)
- **GPU**: CUDA 11.0+ (optional, for acceleration)

## 🔧 Installation Steps

### 1. Install RXNRECer

```bash
pip install rxnrecer
```

### 2. Download Data and Models

```bash
# Download all required files (~35.8GB total)
rxnrecer-download-data

# Or download separately
rxnrecer-download-data --data-only      # ~8.8GB
rxnrecer-download-data --models-only    # ~14GB
rxnrecer-download-data --extools-only   # ~13GB
```

### 3. Set Environment Variables (Optional, for S3 mode)

```bash
export LLM_API_KEY="your_api_key_here"
export LLM_API_URL="https://openrouter.ai/api/v1"
```

## 🚀 Usage

### Command Line Interface

```bash
# Basic prediction (S1 mode)
rxnrecer -i input.fasta -o output.tsv -m s1

# Enhanced prediction (S2 mode)
rxnrecer -i input.fasta -o output.tsv -m s2

# Advanced prediction (S3 mode)
rxnrecer -i input.fasta -o output.json -m s3 -f json

# Get help
rxnrecer --help
```

### Python API

```python
import rxnrecer

# Use the main prediction function
success = rxnrecer.predict(
    input_file="input.fasta",
    output_file="output.tsv",
    mode="s1",
    format="tsv",
    batch_size=100
)
```

### Jupyter Notebook

```python
from rxnrecer.config import config as cfg

# Set your API credentials
cfg.LLM_API_KEY = "your_api_key_here"
cfg.LLM_API_URL = "your_api_url_here"
```

## 🔍 Verification

```bash
# Check installation
rxnrecer --version

# Test with sample data
rxnrecer -i data/sample/sample10.fasta -o test_output.tsv -m s1
```

## 🆕 Features

### Smart Caching

```bash
# Check cache status
rxnrecer-cache status

# Clear cache
rxnrecer-cache clear --all
```

### Data Management

```bash
# Download data
rxnrecer-download-data

# Force re-download
rxnrecer-download-data --force
```

## 🔄 Updating

```bash
pip install --upgrade rxnrecer
```

## 🐛 Troubleshooting

### Common Issues

- **Import Errors**: `pip uninstall rxnrecer && pip install rxnrecer`
- **Missing Dependencies**: Dependencies are automatically installed
- **CUDA Issues**: Install CPU-only PyTorch if needed

### Get Help

- **Documentation**: https://github.com/kingstdio/RXNRECer#readme
- **Contact**: zhenkun.shi@tib.cas.cn

---

**Get started: `pip install rxnrecer`**

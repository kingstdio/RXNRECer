# RXNRECer Installation Guide

This guide explains how to install and use RXNRECer using pip.

## 🚀 Quick Installation

### Option 1: Install from PyPI (Recommended)

```bash
# Install the latest stable version
pip install rxnrecer

# Install with development dependencies
pip install rxnrecer[dev]

# Install with all dependencies
pip install rxnrecer[full]
```

### Option 2: Install from GitHub

```bash
# Install directly from GitHub
pip install git+https://github.com/kingstdio/RXNRECer.git

# Install specific branch
pip install git+https://github.com/kingstdio/RXNRECer.git@release

# Install specific version
pip install git+https://github.com/kingstdio/RXNRECer.git@v1.0.0
```

### Option 3: Install in Development Mode

```bash
# Clone the repository
git clone https://github.com/kingstdio/RXNRECer.git
cd RXNRECer

# Install in development mode
pip install -e .

# Install with development dependencies
pip install -e .[dev]
```

## 📋 Prerequisites

### System Requirements
- **Python**: 3.10 or higher
- **RAM**: 32GB+ (recommended for large datasets)
- **Disk Space**: 40GB+ (for models and data)
- **GPU**: CUDA 11.0+ (optional, for acceleration)

### Python Dependencies
RXNRECer will automatically install these dependencies:
- PyTorch 2.0+
- NumPy, Pandas, Scikit-learn
- Transformers, Fair-ESM
- Biopython, RDKit
- Click, Tqdm, PyYAML

## 🔧 Installation Steps

### 1. Create Virtual Environment (Recommended)

```bash
# Using conda
conda create -n rxnrecer python=3.10
conda activate rxnrecer

# Using venv
python -m venv rxnrecer-env
source rxnrecer-env/bin/activate  # Linux/Mac
# or
rxnrecer-env\Scripts\activate     # Windows
```

### 2. Install RXNRECer

```bash
# Install from PyPI
pip install rxnrecer

# Verify installation
python -c "import rxnrecer; print(rxnrecer.__version__)"
```

### 3. Download Data and Models

```bash
# Download data files (~8.6GB)
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz"
tar -xzf data.tar.gz
rm data.tar.gz

# Download model files (~11.9GB)
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz"
tar -xzf ckpt.tar.gz
rm ckpt.tar.gz
```

### 4. Set Environment Variables

```bash
# Set LLM API key for S3 mode (optional)
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

# Access configuration
from rxnrecer.config import config
print(f"Data directory: {config.DATA_DIR}")
```

### Jupyter Notebook

```python
# In a Jupyter notebook
from rxnrecer.config import config
config.LLM_API_KEY = "your_api_key_here"
config.LLM_API_URL = "https://openrouter.ai/api/v1"

# Run prediction
import rxnrecer
rxnrecer.predict("input.fasta", "output.tsv", mode="s1")
```

## 📁 Project Structure After Installation

```
your_project/
├── input.fasta              # Your input protein sequences
├── output.tsv               # Prediction results
├── data/                    # Downloaded data files
│   ├── sample/              # Sample data
│   ├── uniprot/             # UniProt data
│   ├── rhea/                # Rhea database
│   └── ...
└── ckpt/                    # Downloaded model files
    ├── rxnrecer/            # RXNRECer models
    └── prostt5/             # ProSTT5 models
```

## 🔍 Verification

### Check Installation

```bash
# Check if rxnrecer command is available
which rxnrecer
rxnrecer --version

# Check Python import
python -c "import rxnrecer; print('RXNRECer imported successfully')"
```

### Test with Sample Data

```bash
# Run a quick test
rxnrecer -i data/sample/sample10.fasta -o test_output.tsv -m s1

# Check output
head test_output.tsv
```

## 🐛 Troubleshooting

### Common Issues

#### 1. Import Errors
```bash
# Reinstall the package
pip uninstall rxnrecer
pip install rxnrecer
```

#### 2. Missing Dependencies
```bash
# Install missing dependencies
pip install -r requirements.txt
```

#### 3. CUDA Issues
```bash
# Install CPU-only PyTorch
pip uninstall torch torchvision
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
```

#### 4. Permission Issues
```bash
# Use user installation
pip install --user rxnrecer
```

### Get Help

- **Documentation**: https://github.com/kingstdio/RXNRECer#readme
- **Issues**: https://github.com/kingstdio/RXNRECer/issues
- **Contact**: zhenkun.shi@tib.cas.cn

## 🔄 Updating

### Update from PyPI

```bash
pip install --upgrade rxnrecer
```

### Update from GitHub

```bash
pip install --upgrade git+https://github.com/kingstdio/RXNRECer.git
```

## 📊 Performance Tips

1. **Use GPU**: Install CUDA-enabled PyTorch for faster inference
2. **Batch Processing**: Use appropriate batch sizes for your hardware
3. **Memory Management**: Close other applications to free up RAM
4. **SSD Storage**: Use SSD for faster data loading

## 🎯 Next Steps

After successful installation:

1. **Download data and models** from AWS S3
2. **Test with sample data** to verify functionality
3. **Configure LLM API** for S3 mode (optional)
4. **Run predictions** on your protein sequences
5. **Explore the documentation** for advanced features

---

**Happy predicting! 🚀**

## 🆕 New Features

### 📥 Automatic Data Download
After installation, you can automatically download required data and model files:

```bash
# Download all files (~20.5GB total)
rxnrecer-download-data

# Download only data files (~8.6GB)
rxnrecer-download-data --data-only

# Download only model files (~11.9GB)
rxnrecer-download-data --models-only

# Force re-download (overwrite existing files)
rxnrecer-download-data --force
```

### 💾 Smart Caching System
RXNRECer includes an intelligent caching system that automatically stores prediction results:

```bash
# Use caching (default behavior)
rxnrecer -i input.fasta -o output.tsv -m s1

# Disable caching for this run
rxnrecer -i input.fasta -o output.tsv -m s1 --no-cache

# Check cache status
rxnrecer-cache status

# View cache information
rxnrecer-cache info

# Clear all cache
rxnrecer-cache clear --all

# Clear old cache files (older than 7 days)
rxnrecer-cache clear --older-than 7
```

**Cache Benefits:**
- ⚡ **Instant repeated predictions** - Identical inputs return cached results immediately
- 💰 **Resource efficiency** - No redundant model inference
- 🔄 **Transparent operation** - Works automatically without user intervention
- 📊 **Smart identification** - MD5 hash ensures parameter-specific caching

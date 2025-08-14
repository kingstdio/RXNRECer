# RXNRECer

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**RXNRECer** is a deep learning-based tool for predicting enzyme-catalyzed reactions from protein sequences using PLM and LLM.

## 🚀 Features

- **Protein Sequence Analysis**: Process protein sequences in FASTA format
- **Deep Learning Models**: ESM-2 embeddings with BGRU and attention
- **Reaction Prediction**: Predict enzyme-catalyzed reactions with confidence scores
- **Multiple Output Formats**: Support for TSV, CSV, and JSON outputs
- **Batch Processing**: Efficient batch processing for large datasets
- **Ensemble Methods**: Support for ensemble predictions
- **Easy-to-use CLI**: Simple command-line interface

## 📋 Requirements

- Python 3.8+
- PyTorch 2.0+
- CUDA (optional, for GPU acceleration)

## 🛠️ Installation

### Quick Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/rxnrecer.git
cd rxnrecer

# Run the setup script
./scripts/setup_environment.sh
```

### Manual Installation

```bash
# Create virtual environment
python -m venv rxnrecer
source rxnrecer/bin/activate  # On Windows: rxnrecer\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

## 🚀 Quick Start

### Basic Usage

```python
import rxnrecer

# Initialize predictor
predictor = rxnrecer.RXNRECerPredictor()

# Load model (if you have a pre-trained model)
predictor.load_model("path/to/model.pth")

# Predict reactions
predictions = predictor.predict("input.fasta")
```

### Command Line Interface

```bash
# Basic prediction
rxnrecer predict input.fasta output.tsv

# With custom parameters
rxnrecer predict input.fasta output.json \
    --batch-size 50 \
    --top-k 10 \
    --format json \
    --equations

# Verbose output
rxnrecer predict input.fasta output.tsv --verbose
```

## 📁 Project Structure

```
rxnrecer/
├── rxnrecer/                 # Main package
│   ├── core/                # Core functionality
│   ├── models/              # Model definitions
│   ├── data/                # Data processing
│   ├── utils/               # Utility functions
│   ├── config/              # Configuration
│   └── cli/                 # Command-line interface
├── scripts/                 # Setup and utility scripts
├── tests/                   # Test suite
├── docs/                    # Documentation
├── data/                    # Data directory
├── results/                 # Results directory
└── notebooks/               # Jupyter notebooks
```

## 🔧 Configuration

### Model Configuration

```python
from rxnrecer.config.model_config import ModelConfig

config = ModelConfig(
    esm_out_dim=1280,
    gru_h_dim=512,
    att_dim=32,
    dropout_rate=0.2,
    batch_size=32,
    top_k=5
)
```

### Environment Variables

```bash
# Set environment variables
export RXNRECER_ENV=production
export LOG_LEVEL=INFO
export DEVICE=cuda
```

## 📊 Input Format

### FASTA File

```
>P12345
MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG
>P67890
MKLIVWALVLAFLACQGAVLGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEKIPVDGIKIVGDMVEV
```

## 📈 Output Format

### TSV Output

```tsv
uniprot_id	top_predictions	top_probabilities	confidence
P12345	RXN-12345;RXN-67890;RXN-11111	0.8500;0.1200;0.0300	0.8500
P67890	RXN-67890;RXN-12345;RXN-22222	0.9200;0.0600;0.0200	0.9200
```

### JSON Output

```json
[
  {
    "uniprot_id": "P12345",
    "top_predictions": ["RXN-12345", "RXN-67890", "RXN-11111"],
    "top_probabilities": [0.8500, 0.1200, 0.0300],
    "confidence": 0.8500
  }
]
```

## 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=rxnrecer

# Run specific test
pytest tests/test_predictor.py
```

## 📚 Documentation

- [API Documentation](docs/api.md)
- [User Guide](docs/user_guide.md)
- [Developer Guide](docs/developer_guide.md)
- [Examples](docs/examples/)

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [ESM-2](https://github.com/facebookresearch/esm) for protein language models
- [PyTorch](https://pytorch.org/) for deep learning framework
- [BioPython](https://biopython.org/) for bioinformatics tools

## 📞 Contact

- **Author**: Zhenkun Shi
- **Email**: zhenkun.shi@tib.cas.cn
- **Project**: [https://github.com/kingstdio/rxnrecer](https://github.com/kingstdio/rxnrecer)

## 🔄 Changelog

### Version 1.0.0
- Complete project restructuring
- Modular architecture with clear separation of concerns
- Improved configuration management
- Enhanced CLI interface
- Comprehensive documentation
- Unit tests and development tools


# RXNRECer

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI](https://img.shields.io/badge/PyPI-rxnrecer-blue.svg)](https://pypi.org/project/rxnrecer/)

**RXNRECer v1.4.1** is a deep learning framework for predicting enzyme-catalyzed reactions directly from protein sequences.
It is the official implementation of  
**“RXNRECer: Active Learning with Protein Language Models for Fine-Grained Enzyme Reaction Prediction.”**

🎉 **Now available on PyPI for easy installation!**

---

## 🚀 Features

- **Multi-stage prediction framework**
  - S1: Direct reaction-level prediction
  - S2: Multi-source ensemble integration
  - S3: LLM-based prompt-guided reasoning and interpretability
- **Protein sequence analysis** from FASTA input
- **Deep learning models** based on ESM-2 protein language models
- **GPU acceleration** with CUDA support
- **Easy-to-use CLI** for large-scale batch inference
- **Smart caching** for efficient repeated predictions
- **Web/API deployment notes** for integrating RXNRECer outputs into a browser or service workflow

---

## 📋 System Requirements

- **Python**: 3.10+
- **NumPy**: 1.21 to <2.0 (current compatibility range)
- **PyTorch/CUDA**: use a CUDA-enabled PyTorch build compatible with the installed NVIDIA GPU and driver; S1 can also run on CPU
- **Memory**: ≥32 GB RAM
- **Disk space**: ≥40 GB (for data and model files)
- **Container runtime**: Apptainer or SingularityCE is required for the complete S2 workflow because the CatFam and ECRECer components are distributed as `.sif` images. It is not required for S1.
- **Stage-specific hardware**: S1, S2, and S3 support both CPU execution and GPU acceleration. For S2 prediction on a GPU, at least 16 GB VRAM is recommended.

### S2 Dependencies

The complete RXNRECer-S2 workflow depends on [Apptainer](https://apptainer.org/) or [SingularityCE](https://sylabs.io/singularity/) to run the bundled CatFam and ECRECer container images:

```text
extools/ec/catfam.sif
extools/ec/ecrecer.sif
```

Download these images with `rxnrecer-download-data --extools-only` and confirm that the `singularity` command is available before running `-m s2`. S1 does not require Apptainer/Singularity. S3 runs S2 before the LLM step and therefore has the same S2 dependencies.

### Compatibility

Dependency compatibility is validated for each release rather than tied to one permanent PyTorch/CUDA combination. Consult the current package metadata and release notes when installing or upgrading, especially on newer GPU architectures.

See the [Installation Guide](docs/INSTALL.md) for the validated environment profile, container-runtime checks, hardware requirements, and reproducibility guidance.

---

## ⏱ Installation Time

- **Typical installation time**: ~10–20 minutes (excluding data/model downloads)
- **Data & model download time**: ~15–30 minutes (network dependent)

---

## 🚀 Quick Start

### 1. Install (Recommended)

```bash
# Install from PyPI
python -m pip install "numpy>=1.21,<2.0" rxnrecer

# Or install the latest version from GitHub
python -m pip install "numpy>=1.21,<2.0" git+https://github.com/kingstdio/RXNRECer.git
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

### S2/S3 Container Setup

The complete S2 workflow runs CatFam and ECRECer from the downloaded `.sif` images. Install [Apptainer](https://apptainer.org/docs/admin/main/installation.html) or [SingularityCE](https://docs.sylabs.io/guides/master/admin-guide/installation.html), then download and check the external tools:

```bash
rxnrecer-download-data --extools-only
singularity --version
singularity inspect extools/ec/catfam.sif
singularity inspect extools/ec/ecrecer.sif
singularity exec --nv extools/ec/ecrecer.sif nvidia-smi
```

S1 does not use these containers. S3 runs S2 before the LLM step and therefore uses the same container setup. See [Container Runtime Setup](docs/INSTALL.md#container-runtime-setup-for-s2s3) for Apptainer compatibility, GPU checks, and troubleshooting.

### 3. Run Prediction

```bash
# Basic prediction (S1 mode)
rxnrecer -i input.fasta -o output.tsv -m s1

# Detailed prediction (S2 mode)
rxnrecer -i input.fasta -o output.tsv -m s2

# LLM reasoning (S3 mode, requires API key)
rxnrecer -i input.fasta -o output.json -m s3 -f json
```

### ⏱ Expected Runtime (Demo)
- S1 / S2 inference: ~1–3 minutes for ~100 proteins (tested on a typical workstation: 32 GB RAM + 1×GPU)
- S3 (LLM reasoning): +5–30 seconds / protein (depends on API latency and chosen LLM)
- Training: not included in the demo
- CPU-only inference: supported but significantly slower (not recommended for large batches)

Note: The demo focuses on inference and usage examples. Full benchmarking from the paper requires extra datasets and scripts described in the Methods section.

## 🌐 Public Webserver

RXNRECer is also available as a public webserver:

**https://rxnrecer.biodesign.ac.cn/**

The webserver provides a browser-based interface for users who want to test RXNRECer without installing the package locally. It supports the same typical workflow: submit protein sequences, run RXNRECer prediction, and review reaction-level outputs and supporting information through the web interface.

For large-scale batch inference, reproducible local runs, or integration into custom pipelines, use the PyPI package through the command-line interface or Python API. See **[Webserver and Local Usage](docs/WEB_API_DEPLOYMENT.md)** for how the public webserver relates to the CLI/Python workflow.

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

For S3 mode, configure an OpenAI-compatible chat-completions endpoint. This can be a hosted provider, an institutional gateway, or a local server such as vLLM or SGLang.

```bash
export LLM_API_KEY="your_api_key_here"
export LLM_API_URL="https://api.example.com/v1"
export LLM_MODEL="provider/model-version"
```

The model name should match the exact identifier exposed by your provider or local server. If `LLM_MODEL` is not set, RXNRECer uses its package default.

For local OpenAI-compatible servers:

```bash
# vLLM or SGLang local server
export LLM_API_URL="http://127.0.0.1:8000/v1"
export LLM_API_KEY="EMPTY"
export LLM_MODEL="local-rxnrecer-llm"
```

See **[LLM API Setup](docs/LLM_API_SETUP.md)** for hosted-provider examples, local vLLM/SGLang setup, model-version selection, reproducibility notes, and troubleshooting. For general provider-side API key and endpoint setup, the AIML API documentation (<https://docs.aimlapi.com/>) is also a useful reference.

S3 explanations are intended as human-readable evidence summaries for expert interpretation. They should not be treated as experimental validation or as ground-truth benchmark labels.

### Jupyter Notebook Setup

```python
import os
from rxnrecer.config import config as cfg

# Set your API credentials
cfg.LLM_API_KEY = "your_api_key_here"
cfg.LLM_API_URL = "https://api.example.com/v1"
cfg.LLM_MODEL = "provider/model-version"
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
- **[LLM API Setup](docs/LLM_API_SETUP.md)** - S3 provider, local vLLM/SGLang, model-version, and reproducibility settings
- **[Webserver and Local Usage](docs/WEB_API_DEPLOYMENT.md)** - Public webserver access and the CLI/Python workflow
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

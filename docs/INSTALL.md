# RXNRECer Installation Guide

**Version 1.4.2**

## 🚀 Quick Installation

### Install from PyPI (Recommended)

```bash
python -m pip install "numpy>=1.21,<2.0" rxnrecer
```

### Install from GitHub

```bash
python -m pip install "numpy>=1.21,<2.0" git+https://github.com/kingstdio/RXNRECer.git
```

## 📋 Prerequisites

- **Python**: 3.10+
- **NumPy**: 1.21 to <2.0
- **RAM**: 32GB+ (recommended)
- **Disk Space**: 40GB+ (for models and data)
- **S1 GPU**: optional; CPU inference is supported but slower
- **Complete S2/S3 workflow**: an NVIDIA GPU with at least 16GB VRAM is recommended
- **Container runtime for S2/S3**: Apptainer or SingularityCE with `singularity` command compatibility

S1 does not require a container runtime. The complete S2 workflow uses the CatFam and ECRECer `.sif` images included in the external-tools download; S3 invokes S2 before requesting the LLM explanation and therefore has the same local requirements.

### Validated Environment Profile

For RXNRECer v1.4.0, the following environment provides a known-good reference configuration:

| Component | Validated version |
| --- | --- |
| Python | 3.10 |
| NumPy | 1.23.5 |
| PyTorch | 2.5.1 |
| CUDA build | 12.1 |
| pandas | 1.5.3 |
| scikit-learn | 1.5.1 |
| Apptainer | 1.5.0 |

These versions document the validated environment for this release; they are not permanent requirements for future releases. A reproducible installation should use:

- a Python version supported by the current package metadata;
- NumPy 1.x for the current release;
- a PyTorch build selected for the installed NVIDIA GPU architecture and driver;
- the dependencies resolved from the same tagged RXNRECer release; and
- Apptainer or SingularityCE when running the complete S2/S3 workflow.

For newer GPU architectures, select a PyTorch/CUDA build that supports the installed hardware and driver, then verify RXNRECer in that environment. This keeps future releases able to follow current hardware support while preserving a concrete reference for the present release. When reproducing or reporting a run, record the resolved environment:

```bash
python -m pip freeze > rxnrecer-environment.txt
python -c "import torch; print(torch.__version__, torch.version.cuda)"
nvidia-smi
singularity --version
```

The package metadata and release notes should be treated as the source of truth for the release being installed. NumPy 2.x is not included in the current supported profile because binary scientific dependencies and external components have been validated against NumPy 1.x.

### Hardware Guidance

The minimum local hardware requirements are:

| Stage | CPU/GPU requirement | Common requirements |
| --- | --- | --- |
| S1 | CPU or GPU | 32GB RAM and 40GB disk space |
| S2 | CPU or GPU; at least 16GB VRAM recommended for GPU execution | 32GB RAM and 40GB disk space |
| S3 | CPU or GPU; same local requirements as S2 | 32GB RAM and 40GB disk space |

A locally deployed LLM requires additional CPU/GPU memory and disk space according to the selected model and inference engine.

### Container Runtime Setup for S2/S3

#### 1. Install a SIF-compatible runtime

Install one of the following on the Linux host that will run RXNRECer:

- [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html)
- [SingularityCE installation guide](https://docs.sylabs.io/guides/master/admin-guide/installation.html)

On a managed workstation or HPC cluster, the runtime may already be installed as a system package or environment module. Check both command names:

```bash
command -v singularity
command -v apptainer
singularity --version
apptainer --version
```

RXNRECer v1.4.0 launches the S2 components through the `singularity` command. Many Apptainer installations provide this compatibility command. If `apptainer --version` succeeds but `singularity --version` does not, ask the system administrator to enable Singularity command compatibility or install SingularityCE before running S2/S3.

#### 2. Download the S2 container images

The external-tools archive contains the CatFam and ECRECer images used by S2:

```bash
rxnrecer-download-data --extools-only
ls -lh extools/ec/catfam.sif extools/ec/ecrecer.sif
```

The images are large, so confirm that both downloads completed before prediction. They should remain under `extools/ec/` relative to the RXNRECer project root.

#### 3. Inspect the images and GPU passthrough

```bash
singularity inspect extools/ec/catfam.sif
singularity inspect extools/ec/ecrecer.sif
singularity exec --nv extools/ec/ecrecer.sif nvidia-smi
```

The final command verifies that the GPU-enabled ECRECer container can access the host NVIDIA driver. CatFam does not require `--nv`; ECRECer does.

#### 4. Run a small S2 check

After downloading the data and model archives, run one small input before starting a larger batch:

```bash
rxnrecer -i data/sample/sample10.fasta -o s2_check.tsv -m s2 --batch_size 1
```

The run log should show the MSA, CatFam, ECRECer, and T5 components. Treat a warning that CatFam or ECRECer failed as an incomplete S2 run and check the runtime, image paths, and GPU passthrough before continuing.

## 🔧 Installation Steps

### 1. Install RXNRECer

```bash
python -m pip install "numpy>=1.21,<2.0" rxnrecer
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
export LLM_API_URL="https://api.example.com/v1"
export LLM_MODEL="provider/model-version"
```

For hosted providers, local vLLM/SGLang servers, and reproducible model-version settings, see [LLM API Setup](LLM_API_SETUP.md). For general provider-side API key and endpoint setup, the AIML API documentation (<https://docs.aimlapi.com/>) is also a useful reference.

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
cfg.LLM_API_URL = "https://api.example.com/v1"
cfg.LLM_MODEL = "provider/model-version"
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
- **S2 container command not found**: Install Apptainer or SingularityCE and verify that `singularity --version` succeeds
- **Insufficient GPU memory for S2/S3**: Use an NVIDIA GPU with at least 16GB VRAM or run on CPU.
- **No suitable local GPU**: Use the public webserver at <https://rxnrecer.biodesign.ac.cn/> for interactive testing

### Get Help

- **Documentation**: https://github.com/kingstdio/RXNRECer#readme
- **Contact**: zhenkun.shi@tib.cas.cn

---

**Get started: `pip install rxnrecer`**

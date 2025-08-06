#!/bin/bash

# RXNRECer Installation Script
# This script helps set up the RXNRECer environment

set -e  # Exit on any error

echo "🚀 RXNRECer Installation Script"
echo "================================"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda is not installed. Please install Miniconda or Anaconda first."
    echo "   Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if CUDA is available
if command -v nvidia-smi &> /dev/null; then
    echo "✅ CUDA detected"
    CUDA_AVAILABLE=true
else
    echo "⚠️  CUDA not detected. Will install CPU-only PyTorch."
    CUDA_AVAILABLE=false
fi

# Create conda environment
echo "📦 Creating conda environment 'rxnrecer'..."
conda create -n rxnrecer python=3.10 -y

# Activate environment
echo "🔧 Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rxnrecer

# Install PyTorch based on CUDA availability
if [ "$CUDA_AVAILABLE" = true ]; then
    echo "🔥 Installing PyTorch with CUDA support..."
    conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y
else
    echo "💻 Installing PyTorch (CPU-only)..."
    conda install pytorch torchvision torchaudio cpuonly -c pytorch -y
fi

# Install other dependencies
echo "📚 Installing other dependencies..."
pip install -r requirements.txt

# Install additional conda packages
echo "🔬 Installing bioinformatics packages..."
conda install -c conda-forge rdkit biopython -y

echo ""
echo "✅ Installation completed successfully!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate rxnrecer"
echo ""
echo "To test the installation, run:"
echo "  python test_basic.py"
echo ""
echo "To run a prediction, use:"
echo "  python run_prediction.py <input.fasta> <output.json>" 
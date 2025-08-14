#!/bin/bash

# RXNRECer Environment Setup Script
# This script sets up a complete development environment for RXNRECer

set -e  # Exit on any error

echo "🚀 Setting up RXNRECer development environment..."

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if conda is available
if command -v conda &> /dev/null; then
    print_status "Conda found, using conda environment"
    USE_CONDA=true
else
    print_status "Conda not found, using venv"
    USE_CONDA=false
fi

# Environment name
ENV_NAME="rxnrecer"

# Create virtual environment
if [ "$USE_CONDA" = true ]; then
    print_status "Creating conda environment: $ENV_NAME"
    conda create -n $ENV_NAME python=3.10 -y
    print_success "Conda environment created"
    
    print_status "Activating conda environment"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $ENV_NAME
else
    print_status "Creating virtual environment: $ENV_NAME"
    python3 -m venv $ENV_NAME
    print_success "Virtual environment created"
    
    print_status "Activating virtual environment"
    source $ENV_NAME/bin/activate
fi

# Upgrade pip
print_status "Upgrading pip"
pip install --upgrade pip

# Install PyTorch (with CUDA support if available)
print_status "Installing PyTorch"
if command -v nvidia-smi &> /dev/null; then
    print_status "CUDA detected, installing PyTorch with CUDA support"
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
else
    print_status "CUDA not detected, installing CPU-only PyTorch"
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
fi

# Install core requirements
print_status "Installing core requirements"
pip install -r requirements.txt

# Install development requirements (optional)
read -p "Install development dependencies? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    print_status "Installing development requirements"
    pip install -r requirements-dev.txt
fi

# Install package in development mode
print_status "Installing RXNRECer in development mode"
pip install -e .

# Create necessary directories
print_status "Creating project directories"
mkdir -p data/{raw,processed,models,samples}
mkdir -p results/{predictions,evaluations,logs}
mkdir -p notebooks
mkdir -p tests

# Download sample data (if available)
print_status "Setting up sample data"
if [ -d "data/sample" ]; then
    cp -r data/sample/* data/samples/
    print_success "Sample data copied"
else
    print_warning "Sample data directory not found"
fi

# Create .env file for environment variables
print_status "Creating environment configuration"
cat > .env << EOF
# RXNRECer Environment Configuration
RXNRECER_ENV=development
LOG_LEVEL=INFO
DEVICE=auto

# Data paths
DATA_ROOT=./data
RESULTS_ROOT=./results

# Model configuration
MODEL_CONFIG=default
BATCH_SIZE=32
TOP_K=5
EOF

print_success "Environment configuration created"

# Test installation
print_status "Testing installation"
python -c "import rxnrecer; print('RXNRECer imported successfully')"
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
python -c "import esm; print('ESM imported successfully')"

print_success "Installation test passed!"

# Print next steps
echo
print_success "🎉 RXNRECer environment setup completed!"
echo
echo "Next steps:"
echo "1. Activate the environment:"
if [ "$USE_CONDA" = true ]; then
    echo "   conda activate $ENV_NAME"
else
    echo "   source $ENV_NAME/bin/activate"
fi
echo
echo "2. Test the installation:"
echo "   python -c \"import rxnrecer; print(rxnrecer.get_version())\""
echo
echo "3. Run a quick test:"
echo "   python test_basic.py"
echo
echo "4. Start using RXNRECer:"
echo "   rxnrecer predict input.fasta output.tsv"
echo
echo "For more information, see the documentation in docs/"

print_success "Setup complete! 🚀"

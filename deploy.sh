#!/bin/bash

# RXNRECer Deployment Script
# This script helps you quickly set up and run RXNRECer

set -e

echo "🚀 RXNRECer Deployment Script"
echo "=============================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda is not installed. Please install Miniconda or Anaconda first."
    echo "   Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if git is available
if ! command -v git &> /dev/null; then
    echo "❌ Git is not installed. Please install Git first."
    exit 1
fi

# Create conda environment
echo "📦 Creating conda environment..."
if conda env list | grep -q "rxnrecer"; then
    echo "⚠️  Environment 'rxnrecer' already exists. Updating..."
    conda env update -f environment_rxnrecer-release.yml
else
    echo "🆕 Creating new environment 'rxnrecer'..."
    conda env create -f environment_rxnrecer-release.yml
fi

# Activate environment
echo "🔧 Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rxnrecer

# Install package in development mode
echo "📥 Installing RXNRECer in development mode..."
pip install -e .

# Download data and model files
echo "📥 Downloading data and model files..."
if [ ! -d "data" ] || [ -z "$(ls -A data 2>/dev/null)" ]; then
    echo "   Downloading data files (~8.6GB)..."
    wget -O data.tar.gz "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz"
    tar -xzf data.tar.gz
    rm data.tar.gz
else
    echo "   Data directory already exists, skipping download."
fi

if [ ! -d "ckpt" ] || [ -z "$(ls -A ckpt 2>/dev/null)" ]; then
    echo "   Downloading model files (~11.9GB)..."
    wget -O ckpt.tar.gz "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz"
    tar -xzf ckpt.tar.gz
    rm ckpt.tar.gz
else
    echo "   Model directory already exists, skipping download."
fi

# Test installation
echo "🧪 Testing installation..."
if command -v rxnrecer &> /dev/null; then
    echo "✅ RXNRECer CLI installed successfully!"
    echo "   Version: $(rxnrecer --version 2>/dev/null || echo '1.0.0')"
else
    echo "❌ RXNRECer CLI installation failed."
    exit 1
fi

# Set up environment variables
echo "🔑 Setting up environment variables..."
if [ -z "$LLM_API_KEY" ]; then
    echo "⚠️  LLM_API_KEY not set. You'll need to set it for S3 mode:"
    echo "   export LLM_API_KEY='your_api_key_here'"
fi

# Create sample run script
echo "📝 Creating sample run script..."
cat > run_sample.sh << 'EOF'
#!/bin/bash
# Sample script to run RXNRECer

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rxnrecer

# Set your API key (required for S3 mode)
# export LLM_API_KEY="your_api_key_here"

# Run basic prediction (S1 mode)
echo "Running S1 prediction..."
rxnrecer -i data/sample/sample10.fasta -o results/s1_output.tsv -m s1

# Run enhanced prediction (S2 mode)
echo "Running S2 prediction..."
rxnrecer -i data/sample/sample10.fasta -o results/s2_output.tsv -m s2

# Run advanced prediction (S3 mode) - requires API key
if [ ! -z "$LLM_API_KEY" ]; then
    echo "Running S3 prediction..."
    rxnrecer -i data/sample/sample10.fasta -o results/s3_output.json -m s3 -f json
else
    echo "Skipping S3 prediction (LLM_API_KEY not set)"
fi

echo "✅ All predictions completed!"
echo "📁 Results saved in results/ directory"
EOF

chmod +x run_sample.sh

echo ""
echo "🎉 RXNRECer deployment completed successfully!"
echo ""
echo "📋 Next steps:"
echo "   1. Activate the environment: conda activate rxnrecer"
echo "   2. Set your LLM API key: export LLM_API_KEY='your_key'"
echo "   3. Run sample predictions: ./run_sample.sh"
echo "   4. Check results in the results/ directory"
echo ""
echo "📚 For more information, see:"
echo "   - README.md: Installation and usage guide"
echo "   - RELEASE_NOTES.md: Version information"
echo "   - rxnrecer --help: CLI help"
echo ""
echo "🚀 Happy predicting!"

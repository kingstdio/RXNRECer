# Data Directory

This directory contains all the data files required for RXNRECer to function properly.

## 📁 Directory Structure

```
data/
├── cpd_svg/           # Compound SVG representations
├── dict/              # Dictionary files for mapping and lookup
├── datasets/          # Training and validation datasets
├── uniprot/           # UniProt protein sequence data
├── rxn_json/          # Reaction data in JSON format
├── sample/            # Sample data files for testing
├── chebi/             # ChEBI chemical entity data
├── rhea/              # Rhea reaction database
└── feature_bank/      # Pre-computed feature vectors
```

## 📥 Download Instructions

**Important**: Due to GitHub file size limits, you need to download the data files from AWS S3.

### Option 1: Direct Download from AWS S3

```bash
# Download data files (~8.6GB)
aws s3 cp s3://tibd-public-datasets/rxnrecer/data.tar.gz .

# Extract to current directory
tar -xzf data.tar.gz

# Clean up
rm data.tar.gz
```

### Option 2: Using wget/curl

```bash
# Download data files (~8.6GB)
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/data.tar.gz"

# Extract to current directory
tar -xzf data.tar.gz

# Clean up
rm data.tar.gz
```

### Option 3: Using AWS CLI (Recommended)

```bash
# Install AWS CLI if not already installed
pip install awscli

# Configure AWS credentials (if needed)
aws configure

# Download data files
aws s3 cp s3://tibd-public-datasets/rxnrecer/data.tar.gz .

# Extract to current directory
tar -xzf data.tar.gz

# Clean up
rm data.tar.gz
```

## 🔍 Data Contents

### cpd_svg/
- Chemical compound SVG representations
- Used for visualization and structure analysis

### dict/
- Mapping dictionaries for various identifiers
- UniProt to Rhea mappings
- Chemical compound mappings

### datasets/
- Training datasets for model development
- Validation datasets for model evaluation
- Test datasets for final testing

### uniprot/
- UniProt protein sequence data
- Protein annotations and metadata
- Sequence features and domains

### rxn_json/
- Reaction data in JSON format
- Enzyme-catalyzed reaction information
- Substrate and product mappings

### sample/
- Sample data files for testing
- Small datasets for quick validation
- Example input files

### chebi/
- ChEBI chemical entity data
- Chemical compound information
- Molecular structure data

### rhea/
- Rhea reaction database
- Enzyme-catalyzed reaction database
- Reaction classification and hierarchy

### feature_bank/
- Pre-computed feature vectors
- Protein sequence embeddings
- Model input features

## 🚀 Usage

After downloading and extracting the data files:

1. **Verify Installation**:
   ```bash
   ls -la data/
   # Should show all subdirectories with files
   ```

2. **Test with Sample Data**:
   ```bash
   # Use the sample data for testing
   rxnrecer -i data/sample/sample10.fasta -o results/test_output.tsv -m s1
   ```

3. **Check Data Integrity**:
   ```bash
   # Verify data files are accessible
   python -c "from rxnrecer.config import config; print('Data path:', config.DATA_DIR)"
   ```

## 📊 Data Sources

- **UniProt**: Protein sequence and annotation data
- **Rhea**: Enzyme-catalyzed reaction database
- **ChEBI**: Chemical entities of biological interest
- **Custom datasets**: Curated enzyme reaction datasets

## 🔧 Configuration

Data paths are configured in `rxnrecer/config/config.py`:
```python
# Data configuration
DATA_DIR = "data/"
UNIPROT_DIR = "data/uniprot/"
RHEA_DIR = "data/rhea/"
CHEBI_DIR = "data/chebi/"
FEATURE_BANK_DIR = "data/feature_bank/"
```

## ⚠️ Important Notes

- **File Size**: Total data size is approximately 8.6GB
- **Download Time**: Depending on your internet connection, download may take 10-60 minutes
- **Disk Space**: Ensure you have at least 10GB free disk space
- **Permissions**: Make sure you have write permissions in the current directory

## 🐛 Troubleshooting

### Download Issues
- Check your internet connection
- Verify AWS S3 bucket accessibility
- Ensure sufficient disk space

### Extraction Issues
- Verify the downloaded file is complete
- Check file permissions
- Ensure tar command is available

### Permission Issues
- Run `chmod -R 755 data/` after extraction
- Check user permissions on the directory

## 📞 Support

If you encounter issues with data files:
1. Check the download completed successfully
2. Verify file permissions and disk space
3. Check the log files for error messages
4. Contact: zhenkun.shi@tib.cas.cn

## 🔗 Related Links

- **AWS S3 Bucket**: s3://tibd-public-datasets/rxnrecer/
- **UniProt**: https://www.uniprot.org/
- **Rhea**: https://www.rhea-db.org/
- **ChEBI**: https://www.ebi.ac.uk/chebi/

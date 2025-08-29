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

Download data files using the RXNRECer command:

```bash
# Download all data files (~8.6GB)
rxnrecer-download-data --data-only

# Or download everything including models (~20.5GB total)
rxnrecer-download-data
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

After downloading the data files:

1. **Test with Sample Data**:
   ```bash
   rxnrecer -i data/sample/sample10.fasta -o results/test_output.tsv -m s1
   ```

2. **Verify Data Integrity**:
   ```bash
   ls -la data/
   # Should show all subdirectories with files
   ```

## 📊 Data Sources

- **UniProt**: Protein sequence and annotation data
- **Rhea**: Enzyme-catalyzed reaction database
- **ChEBI**: Chemical entities of biological interest
- **Custom datasets**: Curated enzyme reaction datasets

## ⚠️ Important Notes

- **File Size**: Total data size is approximately 8.6GB
- **Download Time**: Depending on your internet connection, download may take 10-60 minutes
- **Disk Space**: Ensure you have at least 10GB free disk space
- **Permissions**: Make sure you have write permissions in the current directory

## 🐛 Troubleshooting

### Download Issues
- Check your internet connection
- Ensure sufficient disk space
- Verify the download completed successfully

### Permission Issues
- Run `chmod -R 755 data/` after extraction
- Check user permissions on the directory


## 🔗 Related Links

- **UniProt**: https://www.uniprot.org/
- **Rhea**: https://www.rhea-db.org/
- **ChEBI**: https://www.ebi.ac.uk/chebi/

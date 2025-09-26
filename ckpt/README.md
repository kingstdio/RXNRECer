# Checkpoints Directory

This directory contains pre-trained model checkpoints for RXNRECer.

## 📁 Directory Structure

- `rxnrecer/` - RXNRECer model checkpoints
- `prostt5/` - ProSTT5 protein language model checkpoints

## 📥 Download Instructions

Download model files using the RXNRECer command:

```bash
# Download all model files (~14GB)
rxnrecer-download-data --models-only

# Or download everything including data (~35.8GB total)
rxnrecer-download-data
```

## 🧠 Model Information

### RXNRECer Models
- **S1 Model**: Basic reaction prediction using ESM-2 embeddings
- **S2 Model**: Enhanced prediction with reaction details  
- **S3 Model**: Advanced prediction with LLM reasoning

### ProSTT5 Models
- **Protein Language Model**: Pre-trained on large protein sequence datasets
- **Embedding Models**: For feature extraction and representation learning

## 🚀 Usage

After downloading, RXNRECer automatically uses these models for:
- Loading pre-trained weights
- Inference and prediction
- Feature extraction

## ⚠️ Important Notes

- Model files are large (~14GB total)
- Ensure sufficient disk space before downloading
- Models are optimized for GPU inference
- CPU inference is supported but slower


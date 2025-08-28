# Checkpoints Directory

This directory contains pre-trained model checkpoints and weights for RXNRECer.

## 📁 Directory Structure

- `rxnrecer/` - RXNRECer model checkpoints
- `prostt5/` - ProSTT5 protein language model checkpoints

## 📥 Download Instructions

Due to GitHub file size limits, you need to download the model files separately:

```bash
# Download model files (~11.9GB)
wget "https://tibd-public-datasets.s3.us-east-1.amazonaws.com/rxnrecer/ckpt.tar.gz"

# Extract to current directory
tar -xzf ckpt.tar.gz

# Clean up
rm ckpt.tar.gz
```

## 🔍 File Types

The following file types are ignored by Git to prevent large files from being committed:
- `*.tar.gz`, `*.zip` - Archive files
- `*.h5`, `*.hdf5` - HDF5 model files
- `*.pkl`, `*.pickle` - Pickle model files
- `*.pth`, `*.pt` - PyTorch model files
- `*.ckpt` - Checkpoint files
- `*.model` - Generic model files
- `*.weights` - Weight files
- `*.onnx` - ONNX model files
- `*.pb` - TensorFlow model files
- `*.tflite` - TensorFlow Lite model files

## 🧠 Model Information

### RXNRECer Models
- **S1 Model**: Basic reaction prediction using ESM-2 embeddings
- **S2 Model**: Enhanced prediction with reaction details
- **S3 Model**: Advanced prediction with LLM reasoning

### ProSTT5 Models
- **Protein Language Model**: Pre-trained on large protein sequence datasets
- **Embedding Models**: For feature extraction and representation learning

## 🚀 Usage

After downloading the model files, RXNRECer will automatically use them for:
- Loading pre-trained weights
- Inference and prediction
- Feature extraction
- Model fine-tuning

## 📊 Model Performance

These models have been trained on:
- Large-scale enzyme reaction datasets
- Diverse protein sequence collections
- Validated on benchmark datasets
- Optimized for production use

## 🔧 Configuration

Model paths are configured in `rxnrecer/config/config.py`:
```python
# Model checkpoint paths
CHECKPOINT_DIR = "ckpt/"
RXNRECER_MODEL_PATH = "ckpt/rxnrecer/"
PROSTT5_MODEL_PATH = "ckpt/prostt5/"
```

## 📞 Support

If you encounter issues with model files:
1. Verify the download completed successfully
2. Check file permissions and disk space
3. Ensure PyTorch version compatibility
4. Contact: zhenkun.shi@tib.cas.cn

## ⚠️ Important Notes

- Model files are large (~11.9GB total)
- Ensure sufficient disk space before downloading
- Models are optimized for GPU inference
- CPU inference is supported but slower

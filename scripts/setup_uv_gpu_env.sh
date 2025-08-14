#!/usr/bin/env bash
set -euo pipefail

echo "[rxnrecer] Setting up project environment with uv (GPU-enabled if available)..."

# Ensure common local bin paths are on PATH (uv installer uses these)
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"

# Enter project root (script is expected to live in scripts/)
PROJECT_ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$PROJECT_ROOT"

VENV_DIR="$PROJECT_ROOT/.venv"
PYTHON_BIN="$VENV_DIR/bin/python"

echo "[rxnrecer] Project root: $PROJECT_ROOT"

# 1) Ensure uv is installed
if ! command -v uv >/dev/null 2>&1; then
  echo "[rxnrecer] uv not found. Installing uv..."
  curl -Ls https://astral.sh/uv/install.sh | sh -s -- --yes
  export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
else
  echo "[rxnrecer] uv found: $(uv --version)"
fi

# 2) Create virtual environment with uv (non-interactive; clear if exists)
echo "[rxnrecer] Creating virtual environment at $VENV_DIR (will clear if exists)"
UV_VENV_CLEAR=1 uv venv "$VENV_DIR"

# 3) Detect CUDA and choose appropriate PyTorch wheel index
WHEEL_TAG="cu118" # sensible default
if command -v nvidia-smi >/dev/null 2>&1; then
  CUDA_MINOR=$(nvidia-smi 2>/dev/null | sed -n 's/.*CUDA Version: \([0-9]\+\.[0-9]\+\).*/\1/p' | head -n1 || true)
  case "${CUDA_MINOR:-}" in
    12.4*) WHEEL_TAG=cu124 ;;
    12.*)  WHEEL_TAG=cu121 ;;
    11.8*) WHEEL_TAG=cu118 ;;
    *)     WHEEL_TAG=cu118 ;;
  esac
  echo "[rxnrecer] Detected CUDA ${CUDA_MINOR:-unknown}; using PyTorch wheel tag: $WHEEL_TAG"
else
  echo "[rxnrecer] nvidia-smi not found; defaulting PyTorch wheel tag to $WHEEL_TAG"
fi

# 4) Install PyTorch matching the CUDA tag
echo "[rxnrecer] Installing PyTorch (tag: $WHEEL_TAG)"
uv pip install --python "$PYTHON_BIN" torch torchvision torchaudio --index-url "https://download.pytorch.org/whl/$WHEEL_TAG"

# 5) Install remaining requirements (exclude torch family and rdkit, install rdkit-pypi)
TMP_REQ=".tmp.requirements.nogpu.txt"
echo "[rxnrecer] Preparing trimmed requirements"
grep -vE '^(torch(|vision|audio)|rdkit)($|[<>= ])' requirements.txt > "$TMP_REQ" || true

echo "[rxnrecer] Installing project Python dependencies via uv"
uv pip install --python "$PYTHON_BIN" -r "$TMP_REQ"
uv pip install --python "$PYTHON_BIN" rdkit-pypi

# 6) Editable install of the project
echo "[rxnrecer] Installing rxnrecer in editable mode"
uv pip install --python "$PYTHON_BIN" -e "$PROJECT_ROOT"

# 7) Jupyter stack (optional but recommended for notebooks)
echo "[rxnrecer] Ensuring Jupyter components are available in this venv"
uv pip install --python "$PYTHON_BIN" -U jupyter-core jupyter-client ipykernel
"$PYTHON_BIN" -m ipykernel install --user --name venv-rxnrecer-release --display-name "venv-rxnrecer-release" || true

# 8) Smoke test (versions and CUDA availability)
echo "[rxnrecer] Running smoke test"
"$PYTHON_BIN" - <<'PY'
import torch, rxnrecer
print("rxnrecer", rxnrecer.get_version())
print("torch", torch.__version__)
print("cuda_available", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device_count", torch.cuda.device_count())
    print("device_name_0", torch.cuda.get_device_name(0))
PY

echo "[rxnrecer] Environment setup complete. You can activate it with:"
echo "  source $VENV_DIR/bin/activate"

# Optional quick run of example script (comment out if not needed)
echo "[rxnrecer] Running examples/basic_usage.py (CPU-bound example)"
"$PYTHON_BIN" examples/basic_usage.py || true

echo "[rxnrecer] Done."



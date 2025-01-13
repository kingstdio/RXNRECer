#!/bin/bash
#SBATCH --job-name=rnxrecer-smi-10fold
#SBATCH --partition=qcpu_23i
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --output=slurmlog/RXNRECer-SMI-10FOLD-%A_%a.out
#SBATCH --error=slurmlog/RXNRECer-SMI-10FOLD-%A_%a.err
#SBATCH --array=1-10

# Define embedding method
EMBD_METHOD="tdit5"  # Using a single embedding method for all folds

# Fold number is directly based on SLURM_ARRAY_TASK_ID
FOLD_NUM=$SLURM_ARRAY_TASK_ID

# Output task start information
echo "Job started: $(date)"
echo "Running for fold_num=${FOLD_NUM} with embedding method=${EMBD_METHOD}"
echo "Task is running on node: $(hostname)"

# Set PYTHONPATH to include the project root
export PYTHONPATH=$PYTHONPATH:/hpcfs/fhome/shizhenkun/codebase/RXNRECer

# Run the Python script with the specified fold number and embedding method
/hpcfs/fhome/shizhenkun/miniconda3/envs/rxnrecer/bin/python \
/hpcfs/fhome/shizhenkun/codebase/RXNRECer/methods/CaclSimiProteinWithfoldNum.py \
--foldnum $FOLD_NUM --embdmethod $EMBD_METHOD

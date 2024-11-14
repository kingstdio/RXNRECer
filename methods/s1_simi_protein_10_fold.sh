#!/bin/bash
#SBATCH --job-name=rnxrecer-smi-10fold
#SBATCH --partition=qcpu_23i
#SBATCH --cpus-per-task=20
#SBATCH --mem=150G
#SBATCH --output=RXNRECer-SMI-10FOLD-%A_%a.out
#SBATCH --error=RXNRECer-SMI-10FOLD-%A_%a.err
#SBATCH --array=0-29

# Define fold numbers and embedding methods
FOLD_NUMS=(1 2 3 4 5 6 7 8 9 10)
EMBD_METHODS=("esm" "unirep" "t5")

# Calculate fold and embedding method based on the task ID
FOLD_IDX=$(( SLURM_ARRAY_TASK_ID / 3 ))
EMBD_IDX=$(( SLURM_ARRAY_TASK_ID % 3 ))
FOLD_NUM=${FOLD_NUMS[$FOLD_IDX]}
EMBD_METHOD=${EMBD_METHODS[$EMBD_IDX]}

# Load the necessary modules (if any)
# Set PYTHONPATH to include the project root
export PYTHONPATH=$PYTHONPATH:/hpcfs/fhome/shizhenkun/codebase/RXNRECer

# Run the Python script with the specified fold number and embedding method
/hpcfs/fhome/shizhenkun/miniconda3/envs/rxnrecer/bin/python /hpcfs/fhome/shizhenkun/codebase/RXNRECer/methods/simi_protein_10_fold.py --foldnum $FOLD_NUM --embdmethod $EMBD_METHOD
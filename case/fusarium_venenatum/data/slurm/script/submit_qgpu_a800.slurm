#!/bin/bash
#SBATCH --job-name=esmfold-qgpu_a800
#SBATCH --partition=qgpu_a800
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --gres=gpu:1
#SBATCH --output=/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/fusarium_venenatum/data/slurm/log/esmfold-qgpu_a800-%A_%a.out
#SBATCH --error=/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/fusarium_venenatum/data/slurm/log/esmfold-qgpu_a800-%A_%a.err
#SBATCH --array=1-12

INPUT_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/case/fusarium_venenatum/data/slurm/temp"
OUTPUT_DIR="/hpcfs/fpublic/database/pdbs/esmfold/pdbs"
INPUT_FILE="qgpu_a800_chunk_${SLURM_ARRAY_TASK_ID}.fasta"
EXTRA_ARGS=""
case $SLURM_ARRAY_TASK_ID in
    1) EXTRA_ARGS='--chunk-size=128' ;;
    2) EXTRA_ARGS='--chunk-size=128' ;;
    3) EXTRA_ARGS='--chunk-size=128' ;;
    4) EXTRA_ARGS='--chunk-size=128' ;;
    5) EXTRA_ARGS='--chunk-size=128' ;;
    6) EXTRA_ARGS='--chunk-size=128' ;;
    7) EXTRA_ARGS='--chunk-size=128' ;;
    8) EXTRA_ARGS='--chunk-size=128' ;;
    9) EXTRA_ARGS='--chunk-size=128' ;;
    10) EXTRA_ARGS='--chunk-size=128' ;;
    11) EXTRA_ARGS='--chunk-size=128' ;;
    12) EXTRA_ARGS='--chunk-size=64' ;;
esac

mkdir -p $OUTPUT_DIR

echo "[`date`] Running on node: `hostname`"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running command:"
echo "singularity exec --nv --bind /hpcfs:/hpcfs /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \
    bash /esmfold.sh -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS"

singularity exec --nv --bind /hpcfs:/hpcfs /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \
    bash /esmfold.sh -i $INPUT_DIR/$INPUT_FILE -o $OUTPUT_DIR $EXTRA_ARGS

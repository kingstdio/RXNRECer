#!/bin/bash
#SBATCH --job-name=embed-prostT5-3di    # Job name
#SBATCH --partition=qgpu_3090           # Partition name
#SBATCH --cpus-per-task=8              # Number of CPUs per task
#SBATCH --mem=64G                       # Memory allocation
#SBATCH --gres=gpu:1                    # Number of GPUs per task
#SBATCH --output=slurmlog/embed-3di-%A_%a.out  # Standard output log
#SBATCH --error=slurmlog/embed-3di-%A_%a.err   # Standard error log
#SBATCH --array=1-32                    # Task array index (1 to 16)

# Dynamic input/output directories
INPUT_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/intermediate/foldseek/3diembd"
OUTPUT_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/intermediate/foldseek/3diembd"
INPUT_FILE="3difold${SLURM_ARRAY_TASK_ID}.fasta"
OUTPUT_FILE="${OUTPUT_DIR}/3difold${SLURM_ARRAY_TASK_ID}.feather"

# Make sure the output directory exists
mkdir -p $OUTPUT_DIR

# Print job info
echo "Processing input file: ${INPUT_FILE}"
echo "Saving embeddings to: ${OUTPUT_FILE}"

# 输出任务开始信息
echo "Job started: $(date)"
echo "Running foldseek for fold_num=${FOLD_NUM} with ${SLURM_CPUS_PER_TASK} threads"
echo "Task is running on node: $(hostname)"
echo "Temporary directory: ${TEMP_DIR}"


source /etc/profile.d/conda.sh
conda activate /hpcfs/fhome/shizhenkun/miniconda3/envs/rxnrecer

# Execute the Python script
python /hpcfs/fhome/shizhenkun/codebase/RXNRECer/modules/structure/embedProstT5_3di.py \
--input ${INPUT_DIR}/${INPUT_FILE} \
--output ${OUTPUT_FILE} \
--half 1 \
--is_3Di 1 \
--per_protein 1

# Done
echo "Job completed for ${INPUT_FILE}, embeddings saved to ${OUTPUT_FILE}"

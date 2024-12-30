#!/bin/bash
#SBATCH --job-name=rnxrecer-smi-10fold    # 任务名称
#SBATCH --partition=qgpu_a800            # 分区名称
#SBATCH --cpus-per-task=16               # 每个任务使用的 CPU 数量
#SBATCH --mem=120G                        # 每个任务分配的内存
#SBATCH --gres=gpu:1                     # 每个任务分配的 GPU 数量
#SBATCH --output=slurmlog/RXNRECer-SMI-10FOLD-%A_%a.out  # 标准输出文件
#SBATCH --error=slurmlog/RXNRECer-SMI-10FOLD-%A_%a.err   # 标准错误文件
#SBATCH --array=0-9                     # 数组任务索引范围

# 动态设置输入文件和输出目录
INPUT_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/intermediate/esmfold/input"
OUTPUT_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/intermediate/esmfold/output"
INPUT_FILE="chunk${SLURM_ARRAY_TASK_ID}.fasta"

# 确保输出目录存在
mkdir -p $OUTPUT_DIR

# 执行命令
# singularity exec /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \
#     bash /esmfold.sh -i ${INPUT_DIR}/${INPUT_FILE} -o ${OUTPUT_DIR} --cpu-only


singularity exec --nv /hpcfs/fpublic/container/singularity/app/esmfold/esmfold.sif \
bash /esmfold.sh \
-i ${INPUT_DIR}/${INPUT_FILE} \
-o ${OUTPUT_DIR}

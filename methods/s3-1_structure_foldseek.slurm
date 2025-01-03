#!/bin/bash
#SBATCH --job-name=foldseek_batch                       # 作业名称
#SBATCH --partition=qcpu_23i                            # 队列名称
#SBATCH --nodes=1                                       # 每个作业使用一个节点
#SBATCH --ntasks-per-node=1                             # 每个节点的任务数（单任务）
#SBATCH --cpus-per-task=63                              # 每个任务使用的 CPU 核心数
#SBATCH --mem=500G                                      # 每个节点使用的内存
#SBATCH --array=0-9                                     # 数组任务索引范围
#SBATCH --output=slurmlog/RXNRECer-FOLDSEEK-%A_%a.out   # 标准输出文件
#SBATCH --error=slurmlog/RXNRECer-FOLDSEEK-%A_%a.err    # 标准错误文件

# 配置变量
FOLD_NUM=$((SLURM_ARRAY_TASK_ID + 1)) # 将数组索引转换为 fold_num (1 到 10)
CFG_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/intermediate/foldseek"
RESULTS_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/results/structural"
BASE_TEMP_DIR="/hpcfs/fhome/shizhenkun/codebase/RXNRECer/temp"
TEMP_DIR="${BASE_TEMP_DIR}/task_${SLURM_ARRAY_TASK_ID}" # 动态临时目录
FSEEK_COLS="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,prob,lddt,alntmscore,qtmscore,ttmscore"

# 创建日志目录（如果不存在）
mkdir -p slurmlog

# 为当前任务创建独立的临时目录
mkdir -p ${TEMP_DIR}

# 输出任务开始信息
echo "Job started: $(date)"
echo "Running foldseek for fold_num=${FOLD_NUM} with ${SLURM_CPUS_PER_TASK} threads"
echo "Task is running on node: $(hostname)"
echo "Temporary directory: ${TEMP_DIR}"

# 设置 OpenMP 线程数
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# 运行命令
time foldseek easy-search ${CFG_DIR}/fold${FOLD_NUM}/pdb_test/ ${CFG_DIR}/fold${FOLD_NUM}/DB/fold${FOLD_NUM} ${RESULTS_DIR}/foldseek_res_fold${FOLD_NUM}.m8 ${TEMP_DIR} --alignment-type 0 --format-output ${FSEEK_COLS} --tmscore-threshold 0.2 --tmalign-hit-order 0 --max-seqs 1000

# 清理任务临时目录（可选，如果需要保留临时文件，注释掉以下行）
rm -rf ${TEMP_DIR}
echo "Deleted temporary directory: ${TEMP_DIR}"

# 输出任务结束信息
echo "Job ended: $(date)"

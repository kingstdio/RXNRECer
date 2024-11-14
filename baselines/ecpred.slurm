#!/bin/bash
#SBATCH --job-name=ECRECer_ecpred
#SBATCH --partition=qcpu_18i
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --output=/hpcfs/fhome/shizhenkun/codebase/RXNRECer/temp/log/ecpred/ECRECer_ecpred_%A_%a.out
#SBATCH --error=/hpcfs/fhome/shizhenkun/codebase/RXNRECer/temp/log/ecpred/ECRECer_ecpred_%A_%a.err
#SBATCH --array=1-10

# Define fold number based on the task ID
FOLD_NUM=$SLURM_ARRAY_TASK_ID

# Execute ECPred command with the specified fold number
singularity exec /hpcfs/fpublic/container/singularity/app/ecpred/ecpred.sif java -Xmx128G -jar /ECPred/ECPred.jar spmap /hpcfs/fhome/shizhenkun/codebase/RXNRECer/data/datasets/task240524/validation/fold${FOLD_NUM}/valid.fasta /ECPred/ /tmp /hpcfs/fhome/shizhenkun/codebase/RXNRECer/baselines/results/ec_methods/ecpred/res_fold_${FOLD_NUM}.txt

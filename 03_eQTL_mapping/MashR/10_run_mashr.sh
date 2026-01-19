#!/bin/bash
# Job name:
#SBATCH --job-name=mashr
#SBATCH --account=wujie
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-6226
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=1

source Absolute_path.sh
cd ${mashr}
Rscript 3_run_mashr.R \
${mashr}/mashr_input_strong_set.z.txt \
${mashr}/mashr_input_random_set.z.txt \
${mashr}/output \
1000000 \
202511

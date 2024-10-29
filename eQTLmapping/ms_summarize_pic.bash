#!/bin/bash
# Job name:
#SBATCH --job-name=model.selection
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=20G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=2:00:00
#
#SBATCH --qos=scavenger

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0
source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh

mkdir ${model_selection}/picture
mkdir ${model_selection}/summarize 

Rscript ms_summarize_pic.r ${model_selection}/ ${model_selection}/Sus_scrofa.Sscrofa11.1.tss.tes ${model_selection}/PIG_Chrome.txt ${Picture}/ ${model_selection}/summarize/



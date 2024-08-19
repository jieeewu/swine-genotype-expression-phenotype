#!/bin/bash
# Job name:
#SBATCH --job-name=TWAS
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
#SBATCH --mem=60G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=6:00:00
module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript TWASassociation.R ${tissue} ${pheno_file} ${predict_expr} ${output_path} ${test_input_file} ${code_path} ${gene_anno_file} ${eQTL_eGene_summarize} ${suffix}


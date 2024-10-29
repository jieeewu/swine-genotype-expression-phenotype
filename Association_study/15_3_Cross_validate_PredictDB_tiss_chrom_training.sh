#!/bin/bash
# Job name:
#SBATCH --job-name=MR.mash
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
#SBATCH --mem=10G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0

Rscript ${PredictDB}/15_4_gtex_tiss_chrom_training.R ${chr} ${PredictDB_code} ${tissue} ${snp_input} ${gene_input}

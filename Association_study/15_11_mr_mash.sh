#!/bin/bash
# Job name:
#SBATCH --job-name=MR.mash1
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=15
#
# Memory per node:
#SBATCH --mem=100G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=24:00:00
#
#SBATCH --qos=scavenger

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0

source /mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh
Rscript 15_12_mr_mash.r $start $end 15 ${MR_mash}/all_gene/input ${MR_mash}/all_gene/input/cis_genotype ${MR_mash}/all_gene/input/cis_gene_list ${MR_mash}/all_gene/picture ${MR_mash}/all_gene/output


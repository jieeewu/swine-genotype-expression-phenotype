#!/bin/bash
# Job name:
#SBATCH --job-name=PCAforQTL
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per tas k:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=6:00:00
#
#SBATCH --qos=scavenger

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0
Rscript PCAforQTL.R ${tissue} ${gene_tpm_matrix} ${Four_Tissues_group} ${genename} ${Tissue_expression_group} ${tissue}_selected_K.pdf ${PCAForQTL} ${WGS_data}_for_PC_grm_pca200.eigenvec ${tissue}_gene_exp_adj_PCAForQTL.RData ${tissue}_expr_filtered.pdf > ${log}/${tissue}.Rout 2>&1 



















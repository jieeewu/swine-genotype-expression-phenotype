#!/bin/bash
# Job name:
#SBATCH --job-name=association
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Memory per node:
#SBATCH --mem=64G
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=12:00:00

plink --bfile ${WGS_data} --allow-no-sex --assoc --pfilter 0.00001 --pheno ${gene_exp_adj_PCAForQTL} --all-pheno --out ${output} > ${log} 2>&1


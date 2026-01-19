#!/bin/bash
# Job name:
#SBATCH --job-name=genotypePCA
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
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0

#1.keep individuals
plink --noweb --allow-extra-chr --keep-allele-order --bfile "$raw_WGS_data" --keep "$keepindivi_list" --make-bed --out "$WGS_data"

#2. genotype PCA
plink --noweb --allow-extra-chr --keep-allele-order --bfile "$WGS_data" --nonfounders --chr 1-18,23,24 --indep-pairwise 1000 100 0.25 --maf 0.01 --geno 0.2 --make-bed --out "$WGS_data"_for_PC
gcta64 --bfile "$WGS_data"_for_PC --make-grm --out "$WGS_data"_for_PC_grm --threads 10
gcta64 --grm "$WGS_data"_for_PC_grm --pca 200 --out "$WGS_data"_for_PC_grm_pca200 --threads 10

#3. QC
plink --noweb --allow-extra-chr --keep-allele-order --bfile "$WGS_data" --nonfounders --chr 1-18,23,24 --maf 0.01 --geno 0.2 --make-bed --out "$WGS_data"_QC

#4. PCA
Rscript PCA_plot.R "$WGS_data"_for_PC_grm_pca200.eigenval "$WGS_data"_for_PC_grm_pca200.eigenvec ${group} Genotype_PCA_481

#!/bin/bash
# Job name:
#SBATCH --job-name=estimate.empirical.fdr
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=8
#
# Memory per node:
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
#
#SBATCH --qos=scavenger
module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

Rscript eQTLfdr.R ${genename} ${obs.dir} ${perm.dir} 100 0.05 20 ${out.file} > ${log}/eQTLfdr.log 2>&1

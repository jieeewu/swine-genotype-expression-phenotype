#!/bin/bash
# Job name:
#SBATCH --job-name=RNAseq
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
#SBATCH --time=12:00:00
module load hisat2/2.1.0-gcc-4.8.5
for tissue in Mu Li AF BF
do

hisat2 -p 20 --dta --rna-strandness RF -x sus_scrofa_index -1 ${tissue}_500M_clean_1.fq.gz -2 ${tissue}_500M_clean_2.fq.gz 2 > ${log} | samtools sort -@ 20 -O bam -o ${tissue}_500M.align.sorted.bam - >${tissue}_500M.samsort.log 2>&1
done

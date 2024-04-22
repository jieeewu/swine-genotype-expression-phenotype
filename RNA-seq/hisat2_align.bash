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
#SBATCH --mem=20G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=1:00:00

module load hisat2/2.1.0-gcc-4.8.5
module load samtools/1.11-gcc-4.8.5

"$hisat2" -p 20 --dta --rna-strandness RF -x "$sus_scrofa_index" -1 "$fq1" -2 "$fq2" 2 >"$log.hisat2"\
 |"$samtools" sort -@ 10 -O bam -o "$align_output" - > "$log" 2>&1
 
echo "####################################################################"
date
echo ${Sample} align is done
echo "####################################################################"



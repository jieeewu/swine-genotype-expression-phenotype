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
#SBATCH --time=24:00:00

n=0
while read -r sample; do
  n=$((n + 1))
  
  process_sample() {
    zcat "${1}"_1.clean.fq.gz \
    | seqkit sample --threads 10 -p 0.1 -s "$2" \
    | seqkit head --threads 10 -n 1000000 -o "${1}"_1M_1.clean.fq.gz

    zcat "$path1""${1}/${1}"_2.clean.fq.gz \
    | seqkit sample --threads 10 -p 0.1 -s "$2" \
    | seqkit head --threads 10 -n 1000000 -o "${1}"_1M_2.clean.fq.gz
  }
  
  process_sample "$sample" "$n"
done < samplename.txt


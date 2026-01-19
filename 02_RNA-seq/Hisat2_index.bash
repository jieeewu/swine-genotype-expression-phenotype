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

hisat2-build -p 10 --ss "$sus_scrofa.splicesite" --exon "$sus_scrofa.exon" "$Sscrofa11.1_genome" "$sus_scrofa_index" > "$log.hisat2.index"

echo "####################################################################"
date
echo index is done
echo "####################################################################"

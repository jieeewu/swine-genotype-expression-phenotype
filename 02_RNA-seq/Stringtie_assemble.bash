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

###assemble
"$stringtie" -p 24 -e --rf -A "$assemble_output".abund -G "$ref" -o "$assemble_output" -l "$assemble_output" "$align_output" 2>"$log"

echo "####################################################################"
date
echo ${Sample} assemble is done
echo "####################################################################"


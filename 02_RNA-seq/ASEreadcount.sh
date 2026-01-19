#!/bin/bash
# Job name:
#SBATCH --job-name=ASEread
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
#SBATCH --mem=10G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=1:00:00
for geneid in `cat "$path1"LD_genename_tpm200_492.txt`
do
cat "$ref"Sus_scrofa.Sscrofa11.1.gtf|awk '$3=="exon"{print}'|grep "$geneid" >> "$path1"gene_exon_callSNP.gtf
done
cat "$path1"gene_exon_callSNP.gtf|awk '{print$1"\t"$4"\t"$5}'> "$path1"bcftools_range.txt

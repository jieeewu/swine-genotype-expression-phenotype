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

####1.1 extract splice sites from the gtf file
hisat2_extract_splice_sites.py "$ref_path"Sus_scrofa.Sscrofa11.1.104.gtf > "$path2"sus_scrofa.splicesite

####1.2 use hisat2_extract_exons.py to extract the exons from the gtf file.
hisat2_extract_exons.py "$ref_path"Sus_scrofa.Sscrofa11.1.104.gtf > "$path2"sus_scrofa.exon

####1.3 buid index by hisat2_build command using .splicesite and .exon file
hisat2-build -p 24 --snp "$path2"all_sample_called.snp --ss "$path2"sus_scrofa.splicesite --exon "$path2"sus_scrofa.exon "$ref_path"Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa "$path2"sus_scrofa_index > hisat2_index.log 2>&1

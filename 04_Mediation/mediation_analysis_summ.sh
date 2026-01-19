#!/bin/bash
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
#SBATCH --time=4:00:00
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--bfile SNP_clean_IDchanged_QC_for_eQTL \
--extract ${tissue}_snp_uniq.txt \
--keep-allele-order \
--recode vcf-iid \
--out ${mediate_out}/${tissue}_eQTL_info_raw
python3 mediation_analysis.py $tss_file $snp_gene_file $snp_info_file $output_file




#!/bin/bash
# Job name:
#SBATCH --job-name=TWAS
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=20
#
# Memory per node:
#SBATCH --mem=60G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh

PredictDB_code="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/PredictDB-Tutorial-master/code"
snp_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/snp_input"
gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
predictDB_path="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"
Predixcan="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/Predixcan"

for tissue in Mu Li AF BF
do
#awk 'BEGIN{FS="_";OFS="\t"}{print$1,$2}' ${eQTL_egene_summarize}/${tissue}_ms_uniq_cis_eQTL.list > ${eQTL_egene_summarize}/${tissue}_ms_uniq_cis_eQTL_forbcftools.list
#comm -23 <(cat ${WGS_1604}/GWAS_keep_1512indivi.list|sort) <(awk '{print$1}' ${WGS_path}/${tissue}_*_keep_indivi.txt|sort) > ${Predixcan}/${tissue}_test_dataset_indivi_ID.list

/mnt/research/qgg/software/bcftools-1.13/bcftools \
view ${phasing_1501}/1501DLY_SNP_QC_for_QTL_phased.vcf.gz \
-S ${snp_input}/${tissue}_indivi_list_test.txt \
-R ${eQTL_egene_summarize}/${tissue}_ms_uniq_cis_eQTL_forbcftools.list \
-Oz \
-o ${Predixcan}/${tissue}_80indivi_test_dataset_18295SNP_phased_for_Predixcan.vcf.gz \
--threads 20

#--------------------------------recode 012 输出18666SNP的基因型，方便后续查看基因表达与基因型的相关性
# /mnt/research/qgg/software/plink-v1.90b6.18/plink \
# --bfile ${WGS_1604}/991DLY_SNP_clean_IDchanged_QC_for_QTL \
# --extract  ${model_select_path}/all_tissues_uniq_independent_eQTLid.txt \
# --keep-allele-order \
# --recode A \
# --out ${WGS_1604}/991DLY_18666SNP_clean_IDchanged_QC_for_QTL_012

done


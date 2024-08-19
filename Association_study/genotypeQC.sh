#!/bin/bash
# Job name:
#SBATCH --job-name=QC
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=15
#
# Memory per node:
#SBATCH --mem=60G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh

/mnt/research/qgg/software/plink-v1.90b6.18/plink --noweb --allow-extra-chr --keep-allele-order --bfile ${WGS_1604}/1604DLY_SNP_clean_IDchanged --keep ${WGS_1604}/1556_indivi_WGSid.list --make-bed --out ${WGS_1604}/1556DLY_SNP_clean_IDchanged
/mnt/research/qgg/software/plink-v1.90b6.18/plink --noweb --allow-extra-chr --keep-allele-order --bfile ${WGS_1604}/1556DLY_SNP_clean_IDchanged --nonfounders --chr 1-18,23,24 --indep-pairwise 1000 100 0.25 --maf 0.01 --geno 0.2 --make-bed --out ${WGS_1604}/1556DLY_SNP_clean_IDchanged_QC_for_PC
${gcta}/gcta-1.94.1 --bfile ${WGS_1604}/1556DLY_SNP_clean_IDchanged_QC_for_PC --make-grm --out ${WGS_1604}/1556DLY_SNP_clean_IDchanged_QC_for_PC_grm --threads 15
${gcta}/gcta-1.94.1 --grm ${WGS_1604}/1556DLY_SNP_clean_IDchanged_QC_for_PC_grm --pca 200 --out ${WGS_1604}/1556DLY_SNP_clean_IDchanged_QC_for_PC_grm_pca200 --threads 15
/mnt/research/qgg/software/plink-v1.90b6.18/plink --noweb --allow-extra-chr --keep-allele-order --bfile ${WGS_1604}/1556DLY_SNP_clean_IDchanged --keep ${WGS_1604}/keep_1501indivi_rmoutlier_from_genotypePC.list --make-bed --out ${WGS_1604}/1501DLY_SNP_clean_IDchanged
/mnt/research/qgg/software/plink-v1.90b6.18/plink --noweb --allow-extra-chr --keep-allele-order --bfile ${WGS_1604}/1501DLY_SNP_clean_IDchanged --nonfounders --chr 1-18,23,24 --maf 0.01 --geno 0.2 --make-bed --out ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL
/mnt/research/qgg/software/plink-v1.90b6.18/plink --keep-allele-order --bfile ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL --recode vcf-iid --out ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL
/mnt/research/qgg/software/bcftools-1.13/bcftools view -c 0 ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL.vcf -Oz -o ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL.vcf.gz --threads 15 
/mnt/research/qgg/software/bcftools-1.13/bcftools index ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL.vcf.gz
/mnt/research/qgg/software/plink-v1.90b6.18/plink --noweb --allow-extra-chr --keep-allele-order --bfile ${WGS_1604}/1501DLY_SNP_clean_IDchanged --extract ${model_select_path}/all_tissues_uniq_independent_eQTLid.txt --recode vcf-iid --out ${WGS_1604}/1501DLY_18666SNP_clean_IDchanged
/mnt/research/qgg/software/bcftools-1.13/bcftools view -c 0 ${WGS_1604}/1501DLY_18666SNP_clean_IDchanged.vcf -Oz -o ${WGS_1604}/1501DLY_18666SNP_clean_IDchanged.vcf.gz --threads 15
/mnt/research/qgg/software/bcftools-1.13/bcftools index ${WGS_1604}/1501DLY_18666SNP_clean_IDchanged.vcf.gz

/mnt/research/qgg/software/bcftools-1.13/bcftools concat -a -d snp ${WGS_1604}/1501DLY_18666SNP_clean_IDchanged.vcf.gz ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL.vcf.gz -Oz -o ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL_merged.vcf.gz --threads 15
/mnt/research/qgg/software/bcftools-1.13/bcftools index ${WGS_1604}/1501DLY_SNP_clean_IDchanged_QC_for_QTL_merged.vcf.gz

for chr in `seq 1 18`
do
echo -e "pos\tchr\tcM\n 0\t$chr\t0\n 1000000000\t$chr\t1000" > ${phasing_1501}/genMap_1cMperMb_${chr}.txt  #Generate map file
sbatch --export=chr=$chr,WGS_1604=$WGS_1604,phasing_1501=$phasing_1501,log=$log,thread_num=20 phasing.sh
done
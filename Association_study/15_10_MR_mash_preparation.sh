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
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=10G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=48:00:00

source /mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh
cd ${expression_group}
input=${MR_mash}/all_gene/input
cis_genotype=${MR_mash}/all_gene/input/cis_genotype
cis_gene_list=${MR_mash}/all_gene/input/cis_gene_list


#--------------------------------------
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--vcf ${phasing_1501}/1501DLY_SNP_QC_for_QTL_phased.vcf.gz \
--keep-allele-order \
--keep ${WGS_path}/481_indivi_ID_for_plink.list \
--make-bed \
--out ${input}/481_indivi_SNP_QC_for_QTL_phased


#--------------------------------------genotype
#training set
for gene_range in `cat ${expression_group}/all_tissue_9614common_genename.txt|sed '1,450d'` 
do
cat ${input}/Four_tissues_common_gene_cise_1MB_SNP.range |grep ${gene_range} > ${cis_genotype}/${gene_range}_gene_cis.range

/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--bfile ${input}/481_indivi_SNP_QC_for_QTL_phased \
--keep-allele-order \
--extract range ${cis_genotype}/${gene_range}_gene_cis.range \
--maf 0.05 \
--geno 0.05 \
--recode A \
--out ${cis_genotype}/${gene_range}_train_dataset_cisSNP_for_mrmash
done

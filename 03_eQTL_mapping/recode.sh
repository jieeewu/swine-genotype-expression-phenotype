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
#SBATCH --mem=60G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
source Absolute_path.sh

module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

for tissue in Mu Li AF BF
do
/mnt/research/qgg/software/plink-v1.90b6.18/plink \
--vcf ${test_dataset} \
--keep-allele-order \
--recode A \
--out ${test_dataset_recode}
done

sbatch --export=tissue=$tissue,test_snp_input=$test_snp_input,predict_expr_pheno_path=$predict_expr_pheno_path,cor_gene_and_trait_path=$cor_gene_and_trait_path,output=$output,gene_list=$gene_list,expr_file=$expr_file,train_genotype=$train_genotype,test_genotype=$test_genotype,gene_SNP_trio=$gene_SNP_trio TWAS.sh
done

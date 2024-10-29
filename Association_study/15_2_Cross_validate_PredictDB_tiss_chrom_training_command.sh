#!/bin/bash
# Job name:
#SBATCH --job-name=MR.mash
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
#SBATCH --time=4:00:00
PredictDB_code="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/PredictDB-Tutorial-master/code"
snp_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/snp_input"
gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
predictDB_path="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"

for tiss in Mu Li AF BF
do
python ${PredictDB_code}/split_snp_annot_by_chr.py ${snp_input}/${tiss}_annot.txt ${snp_input}/${tiss}_snp_annot
python ${PredictDB_code}/split_genotype_by_chr.py ${snp_input}/${tiss}_genotype_train.txt ${snp_input}/${tiss}_genotype
done

mkdir -p ${predictDB_path}/summary ${predictDB_path}/covariances ${predictDB_path}/weights


for chr in {1..18}
do
sbatch --export=PredictDB_code=$PredictDB_code,chr=$chr,tissue=$tissue,snp_input=$snp_input,gene_input=$gene_input,PredictDB=$predictDB_path ${predictDB_path}/15_3_Cross_validate_PredictDB_tiss_chrom_training.sh
done

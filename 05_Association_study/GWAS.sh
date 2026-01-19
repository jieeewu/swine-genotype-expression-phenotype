#!/bin/bash
# Job name:
#SBATCH --job-name=GWAS.Observe
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
#SBATCH --time=12:00:00

source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh
PC_number=50     
suffix=${PC_number}pc

for breed in S21_3769indivi   
do

mkdir -p ${GWAS_path}
mkdir -p ${GWAS_path}/${breed}
mkdir -p ${GWAS_path}/${breed}/Heritability


${gcta}/gcta-1.94.1 \
--grm ${genotype_PC} \
--pca ${PC_number} \
--out ${genotype_PC}_${suffix} \
--threads 15
cat ${phenotype_path}/${pheno_file} | sed '1d'|awk '{print $1,$2,$3}' > ${GWAS_path}/${breed}/${breed}_gender_cov.txt
cat ${phenotype_path}/${pheno_file}|sed  -n '1p'|tr "\t" "\n"|awk 'NR>3{print}' > ${phenotype_path}/pheno${pheno_count}.list
${gcta}/gcta-1.94.1 \
--bfile ${WGS_path}/${WGS_file} \
--make-grm \
--out ${WGS_path}/${WGS_file} \
--threads 15
${gcta}/gcta-1.94.1 \
--grm ${WGS_path}/${WGS_file} \
--make-bK-sparse 0.05 \
--out ${WGS_path}/${WGS_file}

for pheno in `seq 2 $((pheno_count + 1))`
do
pheno2=`expr $pheno - 1`

pheno_name=`cat ${phenotype_path}/pheno${pheno_count}.list|sed -n ''"$pheno2"'p'`
${gcta}/gcta-1.94.1 \
--bfile ${WGS_path}/${WGS_file} \
--grm-sparse ${WGS_path}/${WGS_file} \
--fastGWA-mlm \
--pheno ${phenotype_path}/${pheno_file} \
--mpheno ${pheno} \
--qcovar ${genotype_PC}_${suffix}.eigenvec \
--covar ${GWAS_path}/${breed}/${breed}_gender_cov.txt \
--thread-num 15 \
--out ${GWAS_path}/${breed}/${breed}_${pheno_name}_${suffix}_gender

done

for pheno in `seq 2 $((pheno_count + 1))`
do
pheno2=`expr $pheno - 1`
pheno_name=`cat ${phenotype_path}/pheno${pheno_count}.list|sed -n ''"$pheno2"'p'`

${gcta}/gcta-1.94.1 \
--reml \
--reml-alg 0 \
--grm ${WGS_path}/${WGS_file} \
--pheno ${phenotype_path}/${pheno_file} \
--mpheno ${pheno} \
--covar ${GWAS_path}/${breed}/${breed}_gender_cov.txt \
--thread-num 15 \
--reml-maxit 100 \
--out ${pheno_name}

echo ${pheno_name} `cat ${pheno_name}.hsq|grep V\(G\)/Vp|cut -d $'\t' -f2,3` `cat ${pheno_name}.hsq|grep Pval | cut -d $'\t' -f2` >>${breed}_pheno${pheno_count}.txt
done
done














#!/bin/bash
# Job name:
#SBATCH --job-name=mashr
#SBATCH --account=wujie
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-6226
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=20

cd ${MASHR}

mkdir -p ./association ./log

ls ./*.snps.txt \
  | sed 's#.*/##;s#.snps.txt##' \
  | parallel -j 20 '
      gene={}

      plink \
	  --bfile '"${WGS}/${tissue}/${WGS_data}"' \
	  --allow-no-sex \
	  --assoc \
	  --pheno '"${pheno_file}"' \
	  --pheno-name ${gene} \
	  --extract ./snp_lists/${gene}.snps.txt \
	  --out ./association/'"${tissue}"'_${gene}_strong \
	  > ./log/'"${tissue}"'_${gene}_strong_log 2>&1'


#!/bin/bash
# Job name:
#SBATCH --job-name=mashr
#SBATCH --account=wujie
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-6226
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=8


cd ${mashr}/random_set
mkdir -p ./association ./log

ls ./snp_lists/*.snps.txt | sed 's#.*/##;s#.snps.txt##'\
  | parallel -j 8 '
      gene={}
      plink \
	  --bfile '"${WGS}/${tissue}/${WGS_file}"' \
	  --allow-no-sex \
	  --assoc \
	  --pheno '"${phenotype}"' \
	  --pheno-name ${gene} \
	  --extract ./snp_lists/${gene}.snps.txt \
	  --out ./association/'"${tissue}"'_${gene}_random \
	  > ./log/'"${tissue}"'_${gene}_random.log 2>&1'


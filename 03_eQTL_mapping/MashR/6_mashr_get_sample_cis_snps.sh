#!/bin/bash
# Job name:
#SBATCH --job-name=mashr
#SBATCH --account=wujie
# Number of nodes needed for use case:
#SBATCH --nodes=1
#SBATCH --partition=cpu-6226
# Tasks per node based on number of cores per node:
#SBATCH --cpus-per-task=1

source Absolute_path.sh
mkdir -p ${mashr}/random_set
cd ${mashr}/random_set

CIS_WINDOW=1000000

awk 'BEGIN{OFS="\t"}{
  print $1, $4-1, $4, $2
}' ${ref_BIM}.bim> ${ref_BIM}.bed

bedtools intersect \
  -a ${OUT_bed} \
  -b ${ref_BIM}.bed \
  -wa -wb > ${OUT_pair}

K=1000
mkdir -p ./snp_lists

gawk -v K=$K '
BEGIN {
  srand(202511)
}
{
  gene = $4
  snp  = $8

  count[gene]++

  if (count[gene] <= K) {
    res[gene, count[gene]] = snp
  } else {
    j = int(rand() * count[gene]) + 1
    if (j <= K) {
      res[gene, j] = snp
    }
  }
}
END {
  for (key in res) {
    split(key, a, SUBSEP)
    gene = a[1]
    print res[key] > ("./snp_lists/" gene ".snps.txt")
  }
}
' ${OUT_pair}


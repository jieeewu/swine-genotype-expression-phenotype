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

cd ${MASHR}

for tissue in Mu Li AF BF
do
sed -n '1p' ${phenotype}|tr " " "\n"|awk 'NR>2' > ${mashr}/${tissue}.list
done
cat ${mashr}/*.list > ${mashr}/all_tested_genes.list

awk 'BEGIN{OFS="\t"}
$3=="gene"{
  match($0,/gene_id "([^"]+)"/,a)
  gene=a[1]

  chr=$1
  strand=$7

  if(strand=="+"){
    tss_start=$4
    tss_end=$4
  } else if(strand=="-"){
    tss_start=$5
    tss_end=$5
  }

  print chr,gene,tss_start,tss_end
}' ${ref}/${gtf_file} > ${mashr}/${TSS}
awk '$1 ~ /^[0-9XY]+$/' ${mashr}/${TSS} > ${mashr}/${TSS_autosomal}

cd ${mashr}

CIS_WINDOW=1000000
BED=${TSS_autosomal}
GENES=./${genes_for_eQTL}.list
OUT_bed=./${genes_for_eQTL}.bed

awk -v w=${CIS_WINDOW} 'BEGIN{OFS="\t"}
NR==FNR {
  g=$1; sub(/\.[0-9]+$/, "", g)
  keep[g]=1
  next
}
{
  g=$2; sub(/\.[0-9]+$/, "", g)
  if (!(g in keep)) next

  start = $3 - w
  if (start < 0) start = 0
  end   = $4 + w

  # BED4: chr start end name
  print $1, start, end, g
}' ${GENES} ${BED} > ${OUT_bed}

#------------------------------------------------------------strong_set
for tissue in BF AF Li
do
python ../1_1_prepare_strong_set.py --tissue ${tissue}
done

python ../1_2_prepare_strong_set.py

wait

awk 'BEGIN{OFS="\t"}
NR>1{
  split($1,a,"_");
  chr=a[1];
  pos=a[2];
  snp=$1;
  gene=$2;
  print chr, pos-1, pos, snp, gene
}' ./strong_pairs.global.txt > ./strong_pairs.snp.bed

bedtools intersect \
  -a ./strong_pairs.snp.bed \
  -b ./all_tested_genes_for_eQTL.bed \
  -wa -wb \
| awk '$5 == $9' \
> ./strong_pairs.cis.bed


mkdir -p ./snp_lists
awk 'BEGIN{OFS="\t"}
{
  gene=$5
  snp=$4
  print snp >> "./snp_lists/"gene".snps.txt"
}
' ./strong_pairs.cis.bed

#----------------------------------------------------------mashr_strong_set_get_association

tissue_list=("Mu" "Li" "AF")
tissue_num_list=("Mu_468" "Li_442" "AF_457")
for i in ${!tissue_list[@]}; do
    tissue=${tissue_list[$i]}
    tissue_num=${tissue_num_list[$i]}
    
    sbatch --export=mashr=$mashr,tissue=$tissue,tissue_num=$tissue_num,WGS=$WGS,pheno_observe=$pheno_observe ../1_3_mashr_strong_set_get_association.sh
done

wait


#----------------------------------------------------------mashr_strong_set_get_input_file
cd ${mashr}

python ../1_4_mashr_strong_set_get_input_file.py



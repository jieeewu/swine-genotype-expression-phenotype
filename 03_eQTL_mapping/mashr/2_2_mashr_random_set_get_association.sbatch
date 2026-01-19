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
cd ${mashr}

tissue_list=("Mu" "Li" "AF" "BF")
tissue_num_list=("Mu_468" "Li_442" "AF_457" "BF_451")
for i in ${!tissue_list[@]}; do
    tissue=${tissue_list[$i]}
    tissue_num=${tissue_num_list[$i]}
    sbatch --export=mashr=$mashr,tissue=$tissue,tissue_num=$tissue_num,WGS=$WGS,pheno_observe=$pheno_observe 2_2_mashr_random_set_get_association.sh
done

wait

python ../2_3_mashr_random_set_get_input.py

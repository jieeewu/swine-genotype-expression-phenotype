#!/bin/bash
# Job name:
#SBATCH --job-name=enrich
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
#SBATCH --mem=20G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

#---------------------------------------------------------------------eQTL_eGene pair
common_gene=/mnt/research/qgg/wu/YFJH/eQTL_mapping/eQTL_output/FDR/common_gene
awk 'NR>1{print $5"_"$6"_"$1}' /Mu_eqtl_cis_trans.txt > ${common_gene}/Muscle_SGEP_eQTL_eGene_pair.list
awk 'NR>1{print $5"_"$6"_"$1}' /Li_eqtl_cis_trans.txt > ${common_gene}/Liver_SGEP_eQTL_eGene_pair.list


awk 'NR>1{print $2"_"$1}' /pigGTEx_Breed_interaction/PigGTEx_v0.significant_eQTL/Muscle.cis_qtl_pairs.significant.txt|cut -d "_" -f1,2,5 > ${common_gene}/Muscle_pigGTEx_eQTL_eGene_pair.list
awk 'NR>1{print $2"_"$1}' /pigGTEx_Breed_interaction/PigGTEx_v0.significant_eQTL/Liver.cis_qtl_pairs.significant.txt|cut -d "_" -f1,2,5 > ${common_gene}/Liver_pigGTEx_eQTL_eGene_pair.list

awk '{print $1}' Ballester/Muscle_D_L_Y.eQTL2 > ${common_gene}/Muscle_Ballester_eQTL_eGene_pair.list
awk '{print $3"_"$5"_"$2}' Ballester/Liver_D_L_Y.eQTL > ${common_gene}/Liver_Ballester_eQTL_eGene_pair.list


module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0

Rscript 12_2_common_gene_and_eQTL_eGene_pair_pic.r
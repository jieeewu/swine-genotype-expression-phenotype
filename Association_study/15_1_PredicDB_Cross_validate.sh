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
#SBATCH --time=12:00:00

PredictDB_code="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/PredictDB-Tutorial-master/code"
snp_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/snp_input"
gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
predictDB_path="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"

for times in 1 2 3 4 5 
do

for tiss in Mu Li AF BF
do
gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
snp_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/snp_input"
cp ${gene_input}/${tiss}_transformed_expression${times}.txt ${gene_input}/${tiss}_transformed_expression.txt
cp ${gene_input}/${tiss}_covariates${times}.txt ${gene_input}/${tiss}_covariates.txt

cp ${snp_input}/${tiss}_genotype_train${times}.txt ${snp_input}/${tiss}_genotype_train.txt
cp ${snp_input}/${tiss}_indivi_list_test${times}.txt ${snp_input}/${tiss}_indivi_list_test.txt
done

sbatch --export=tissue=Mu,times=$times 15_2_Cross_validate_PredictDB_tiss_chrom_training_command.sh
sbatch --export=tissue=Li,times=$times 15_2_Cross_validate_PredictDB_tiss_chrom_training_command.sh
sbatch --export=tissue=AF,times=$times 15_2_Cross_validate_PredictDB_tiss_chrom_training_command.sh
sbatch --export=tissue=BF,times=$times 15_2_Cross_validate_PredictDB_tiss_chrom_training_command.sh

#
jobID=$(sbatch --export=tissue=BF,times=$times 15_2_Cross_validate_PredictDB_tiss_chrom_training_command.sh| cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

#
jobID=$(sbatch --export=times=$times 15_6_Cross_validate_PredictDB_make_db.sh| cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done
#
jobID=$(sbatch --export=times=$times 15_8_keep_individual_and_SNP_prepare_for_Predixcan.sh| cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done
#
jobID=$(sbatch --export=cross=cross$times 15_9_working_Predixcan_for_each_tissue_genoPC_gender.sh| cut -d " " -f 4)
sleep 5m
jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
while [ "$jobState" == "PENDING" ] || [ "$jobState" == "RUNNING" ] || [ "$jobState" == "COMPLETING" ]
do
	sleep 5m
	jobState=$(sacct -j $jobID --format=state | tail -n 1 | tr -d " ")
done

done


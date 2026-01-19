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
#SBATCH --time=2:00:00
module load GCCcore/8.2.0 Python/3.7.2
source /mnt/research/qgg/wu/YFJH/eQTL_mapping/virtualenv_Predixcan/Predixcan/bin/activate
source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh

for tissue in Mu Li AF BF 
do
printf "Predict expression\n\n"
python3 $METAXCAN/Predict.py \
--model_db_path $db/${tissue}_models_filtered_signif.db \
--vcf_genotypes $Predixcan/${tissue}_994indivi_test_dataset_18295SNP_phased_for_Predixcan.vcf.gz \
--vcf_mode genotyped \
--prediction_output ${Predixcan}/Adjust_predict_output/${tissue}_18295SNP_predict.txt \
--prediction_summary_output ${Predixcan}/Adjust_predict_output/${tissue}_18295SNP_summary.txt \
--verbosity 9 \
--throw

printf "association\n\n"

for pheno in `cat ${phenotype}/phenotype.list`
do
python3 $METAXCAN/PrediXcanAssociation.py \
--expression_file ${Predixcan}/Adjust_predict_output/${tissue}_SNP_predict.txt \
--input_phenos_file ${Predixcan}/pheno_for_TWAS.txt \
--input_phenos_column $pheno \
--output ${Predixcan}/Adjust_PrediXcanAssociation/${tissue}/${tissue}_${pheno}_association.txt \
--verbosity 9 \
--throw
done

done
deactivate

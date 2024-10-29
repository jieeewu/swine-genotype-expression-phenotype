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
#SBATCH --mem=20G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=1:00:00
module purge
module load GCC/8.3.0 OpenMPI/3.1.4 R/4.1.0
source /mnt/research/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh

for tissue in Mu Li AF BF
do
mainpath="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"
predictDB_path="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"

mkdir -p ${predictDB_path}/dbs

Rscript ${code}/15_7_PredictDB_make_db.R ${tissue} ${mainpath}
done


#!/bin/bash
# Job name:
#SBATCH --job-name=phasing
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
#SBATCH --time=48:00:00

/mnt/research/qgg/software/shapeit-4.2/bin/shapeit4.2 \
--input ${WGS_file} \
--map ${phasing}/${genMap} \
--region ${chr} \
--output ${phasing_result} \
--thread ${thread_num} \
--log ${log}/phasing/${phasing_log}

/mnt/research/qgg/software/bcftools-1.13/bcftools index ${phasing_result}


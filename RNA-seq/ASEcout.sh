#!/bin/bash
# Job name:
#SBATCH --job-name=ASEread
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
#SBATCH --time=24:00:00
indsfile=$1
start=$2
end=$3

module load bedtools2/2.26.0-gcc-4.8.5
module load samtools/1.11-gcc-4.8.5
module load GATK/4.1.4.1

sample=(`cat "$path3""$indsfile" | awk '{if(NR>='"$start"' && NR<='"$end"')print}'`)
mkdir "$log"${sample[$SLURM_ARRAY_TASK_ID]}
mkdir "$out"${sample[$SLURM_ARRAY_TASK_ID]}
mkdir "$tmp"${sample[$SLURM_ARRAY_TASK_ID]}

bedtools intersect -a "$path1"${sample[$SLURM_ARRAY_TASK_ID]}.align.sorted_22769new.bam -b "$path2"LD_gene_exon_callSNP.gtf > "$bam"${sample[$SLURM_ARRAY_TASK_ID]}.geneexon_callSNP.bam
gatk AddOrReplaceReadGroups\
 --TMP_DIR "$TPM_DIR"\
 -I "$bam"${sample[$SLURM_ARRAY_TASK_ID]}.geneexon_callSNP.bam\
 -O "$tmp"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.addreadgroup.bam\
 -ID 4\
 -LB lib1\
 -PL illumina\
 -PU unit1\
 -SM 20 > "$log"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.add_read_group.log 2>&1

echo $(date) add readgroup ID done 

gatk MarkDuplicates\
 --TMP_DIR "$TPM_DIR"\
 -I "$tmp"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.addreadgroup.bam\
 -O "$out"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.addreadgroup.markdup.bam\
 -M "$metric"${sample[$SLURM_ARRAY_TASK_ID]}.addreadgroup.markdup_metrics.txt\
 --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500\
 --REMOVE_DUPLICATES false\
 --ASSUME_SORTED true\
 --CREATE_INDEX true > "$log"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.markdup.log 2>&1

echo ${sample[$SLURM_ARRAY_TASK_ID]} MarkDuplicates done

mkdir "$log2"${sample[$SLURM_ARRAY_TASK_ID]}
gatk ASEReadCounter\
 --tmp-dir "$TPM_DIR"\
 -R "$ref"Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa\
 -I "$out"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}.addreadgroup.markdup.bam\
 -V "$vcf"IDchanged_bcftools_gene_exon_492.vcf.gz\
 -O "$out2"${sample[$SLURM_ARRAY_TASK_ID]}_ASEreadCount.table > "$log2"${sample[$SLURM_ARRAY_TASK_ID]}/${sample[$SLURM_ARRAY_TASK_ID]}_ASEreadCount.log 2>&1
echo ${sample[$SLURM_ARRAY_TASK_ID]} ASEreadCount done

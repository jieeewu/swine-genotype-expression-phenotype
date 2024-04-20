#!/bin/bash
# Job name:
#SBATCH --job-name=RNAseq
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
#SBATCH --time=24:00:00

"$fastp"/fastp -w 8 -l 150 -q 5 -u 50 -n 15\
 -i "$sample1".fq.gz -I "$sample2".fq.gz\
 -o "$cleandata1".fq.gz -O "$cleandata2".fq.gz\
 --html "$cleandata1".hmtl\
 --json "$cleandata2".json

"$hisat2" -p 20 --dta --rna-strandness RF -x sus_scrofa_index -1 "$fq1" -2 "$fq2" 2 >"$log"\
 |"$samtools" sort -@ 10 -O bam -o "$align_output" - > "$log" 2>&1
 
"$stringtie" -p 6 -e --rf -A "$assemble_output" -G "$ref" -o "$assemble_output" -l "$align_output" 2>"$log"

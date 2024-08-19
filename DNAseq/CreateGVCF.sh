#!/bin/bash
samplefile=$1
start=$2
end=$3

Sample=(`cat $BASH/$samplefile | awk '{if(NR>='"$start"'&&NR<='"$end"')print $1}'`)

## create gvcf
sentieon driver -t 24 -r $REF \
-i ${INPUT}/${Sample[$SLURM_ARRAY_TASK_ID]}_BQSR.bam \
--algo Haplotyper \
--emit_mode GVCF \
${Sample[$SLURM_ARRAY_TASK_ID]}_first.g.vcf.gz



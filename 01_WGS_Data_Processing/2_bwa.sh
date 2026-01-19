#!/bin/bash
source /BIGDATA1/app/toolshs/moduleenv.sh
module load bzip2/1.0.6-gcc-4.8.5
module load bwa/0.7.12-gcc-4.8.5
module load samtools/1.9-gcc-4.8.5

indsfile=$1
start=$2
end=$3

cd $OUTPUT

SampleList=(`cat $BASH/$indsfile | awk '{if(NR>='"$start"'&&NR<='"$end"')print}'`)

RG="@RG\\tID:${SampleList[$SLURM_ARRAY_TASK_ID]}\\tLB:${SampleList[$SLURM_ARRAY_TASK_ID]}_Lib\\tPL:illumina\\tSM:${SampleList[$SLURM_ARRAY_TASK_ID]}"
echo "BAM HEADER IS: ${RG}"

echo "Mapping"
time bwa mem -M -t 24 -R ${RG} $REF ${INPUT}/${SampleList[$SLURM_ARRAY_TASK_ID]}/${SampleList[$SLURM_ARRAY_TASK_ID]}_1.fq.gz ${INPUT}/${SampleList[$SLURM_ARRAY_TASK_ID]}/${SampleList[$SLURM_ARRAY_TASK_ID]}_2.fq.gz | \
samtools view -bSu - | samtools sort -@ 24 -O bam -o ${SampleList[$SLURM_ARRAY_TASK_ID]}_sorted.bam - && echo "** bwa mapping and sorting done **"

samtools index -@ 24 ${SampleList[$SLURM_ARRAY_TASK_ID]}_sorted.bam

echo "####################################################################"
date
echo " ${SampleList[$SLURM_ARRAY_TASK_ID]} is done"
echo "####################################################################"

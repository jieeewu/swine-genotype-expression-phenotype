#!/bin/bash

indsfile=$1
start=$2
end=$3

cd $OUTPUT

Sample=(`cat $BASH/$indsfile | awk '{if(NR>='"$start"'&&NR<='"$end"')print}'`)

## Metric
sentieon driver -t 24 -r $REF \
    -i ${INPUT}/${Sample[$SLURM_ARRAY_TASK_ID]}_sorted.bam \
    --algo GCBias --summary ${Sample[$SLURM_ARRAY_TASK_ID]}.GC_summary.txt ${Sample[$SLURM_ARRAY_TASK_ID]}.GC_metric.txt \
    --algo MeanQualityByCycle ${Sample[$SLURM_ARRAY_TASK_ID]}.MQ_metric.txt \
    --algo QualDistribution ${Sample[$SLURM_ARRAY_TASK_ID]}.QD_metric.txt \
    --algo InsertSizeMetricAlgo ${Sample[$SLURM_ARRAY_TASK_ID]}.IS_metric.txt \
    --algo AlignmentStat ${Sample[$SLURM_ARRAY_TASK_ID]}.ALN_metric.txt
sentieon plot metrics -o ${Sample[$SLURM_ARRAY_TASK_ID]}_metric_plot \
    gc=${Sample[$SLURM_ARRAY_TASK_ID]}.GC_metric.txt \
    mq=${Sample[$SLURM_ARRAY_TASK_ID]}.MQ_metric.txt \
    qd=${Sample[$SLURM_ARRAY_TASK_ID]}.QD_metric.txt \
    isize=${Sample[$SLURM_ARRAY_TASK_ID]}.IS_metric.txt

sentieon driver -t 24 -i ${INPUT}/${Sample[$SLURM_ARRAY_TASK_ID]}_sorted.bam \
    --algo LocusCollector --fun score_info ${Sample[$SLURM_ARRAY_TASK_ID]}.score.gz 
sentieon driver -t 24 -i ${INPUT}/${Sample[$SLURM_ARRAY_TASK_ID]}_sorted.bam \
    --algo Dedup --rmdup --score_info ${Sample[$SLURM_ARRAY_TASK_ID]}.score.gz \
    --metrics ${Sample[$SLURM_ARRAY_TASK_ID]}.dedup_metric.txt ${Sample[$SLURM_ARRAY_TASK_ID]}.deduped.bam

 sentieon driver -t 24 -r $REF -i ${Sample[$SLURM_ARRAY_TASK_ID]}.deduped.bam \
     --algo Realigner \
	 -k /BIGDATA2/scau_jyang_4/index/1.Pig/Ss11.1_dbsnp.nospace.vcf \
     ${Sample[$SLURM_ARRAY_TASK_ID]}_realigned.bam

sentieon driver -t 24 -r $REF \
      -i ${Sample[$SLURM_ARRAY_TASK_ID]}_realigned.bam \
      --algo QualCal \
	  -k /BIGDATA2/scau_jyang_4/index/1.Pig/Ss11.1_dbsnp.nospace.vcf \
	  ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_data.table
sentieon driver -t 24 -r $REF \
      -i ${Sample[$SLURM_ARRAY_TASK_ID]}_realigned.bam \
      -q ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_data.table \
      --algo QualCal \
      -k /BIGDATA2/scau_jyang_4/index/1.Pig/Ss11.1_dbsnp.nospace.vcf ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_data.table.post \
      --algo ReadWriter \
      ${Sample[$SLURM_ARRAY_TASK_ID]}_BQSR.bam
sentieon driver -t 24 \
      --algo QualCal \
      --plot \
      --before ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_data.table \
      --after ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_data.table.post \
      ${Sample[$SLURM_ARRAY_TASK_ID]}_recal.csv 
sentieon plot bqsr -o ${Sample[$SLURM_ARRAY_TASK_ID]}_recal_plots.pdf ${Sample[$SLURM_ARRAY_TASK_ID]}_recal.csv

echo "####################################################################"
date
echo " ${Sample[$SLURM_ARRAY_TASK_ID]} is done"
echo "####################################################################"

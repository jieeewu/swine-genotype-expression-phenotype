#!/bin/bash

samplefile=$1
Sample=(`cat $BASH/$samplefile`)
N=(`cat $samplefile | wc -l`)
((n=$N-1))

sample_gvcf=""
for i in `seq 0 $n`; do
	sample_gvcf=${sample_gvcf}"-v ${INPUT}/${Sample[$i]}_first.g.vcf.gz "
done

sentieon driver -t 48 -r $REF \
--algo GVCFtyper \
-d /BIGDATA2/scau_jyang_4/index/1.Pig/Ss11.1_dbsnp.nospace.vcf \
${sample_gvcf} 498DLY_dbsnp.vcf.gz

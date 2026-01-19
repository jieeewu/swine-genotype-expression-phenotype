#!/bin/bash

source /BIGDATA1/app/toolshs/moduleenv.sh
module load GATK/4.0.2.1
module load bzip2/1.0.6-gcc-4.8.5
module load bcftools/1.3.1-gcc-4.8.5
module load tabix/0.2.6
module load plink/1.90-gcc-4.8.5

OUTPUT=$INPUT

cd $OUTPUT

gatk SelectVariants \
	-V 107licha_dbsnp.vcf.gz \
	-select-type SNP \
	--restrict-alleles-to BIALLELIC \
	-O 107licha_SNP_raw.vcf.gz

gatk VariantFiltration \
        -R ${REF} \
        -V 107licha_SNP_raw.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "filtered" \
        -O 107licha_SNP_mkfilterd.vcf.gz
bcftools view --threads 24 -i 'FILTER="PASS"' 107licha_SNP_mkfilterd.vcf.gz -Oz -o 107licha_SNP_clean.vcf.gz
tabix -p vcf 107licha_SNP_clean.vcf.gz 

plink --vcf 107licha_SNP_cleanSOR.vcf.gz --make-bed --out 107licha_SNP

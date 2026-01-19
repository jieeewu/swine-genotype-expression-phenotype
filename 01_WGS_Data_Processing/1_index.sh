#!/bin/bash
source /BIGDATA1/app/toolshs/moduleenv.sh
module load bwa/0.7.17-gcc-4.8.5
cd /BIGDATA2/scau_jyang_4/index/1.Pig/
bwa index Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

source /BIGDATA1/app/toolshs/moduleenv.sh
module avail samtools
module load samtools/1.11-gcc-4.8.5
cd /BIGDATA2/scau_jyang_4/index/1.Pig/
samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

module load GATK/4.1.4.1
cd /BIGDATA2/scau_jyang_4/index/1.Pig/
gatk CreateSequenceDictionary -R Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

echo "####################################################################"
date
echo " Index is done"
echo "####################################################################"


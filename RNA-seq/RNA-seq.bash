#!/bin/bash

#need to change
main_path="/data3/wujie_data/projects/RNA_seq/zhanglin"
raw="/data3/wujie_data/projects/RNA_seq/temp/wushangkai08017/X101SC21032408-Z01-F001-B1-16/raw"

#usually path
index="/public/home/wujie/RNA_seq/YF_plan/HISAT2/hisat2_index"
ref="/public/home/wujie/RNA_seq/YF_plan/ref"

#software
fastp="/public/home/wujie/software"
hisat2="/public/home/wujie/software/hisat2-2.1.0"
samtools="/public/home/wujie/software/samtools-1.10/bin"
stringtie="/public/home/wujie/software/stringtie-2.1.6.Linux_x86_64"
get_TPM_FPKM="/public/home/wujie/software/stringtie-2.1.6.Linux_x86_64"
#same path
log=${main_path}/log
code=${main_path}/code
cleandata=${main_path}/cleandata
align_output=${main_path}/hisat2
statistic_output=${main_path}/code/statistic_output
assemble_output=${main_path}/stringtie
TPM_FPKM_output=${main_path}/TPM_FPKM


mkdir ${main_path}/cleandata ${main_path}/code ${main_path}/code/statistic_output ${main_path}/hisat2 ${main_path}/log ${main_path}/stringtie ${main_path}/TPM_FPKM

#------------------------------1_fastp.sh
#!/bin/bash
for sample in `cat "$code"/sample.name`
do
mkdir "$cleandata"/"$sample"
"$fastp"/fastp -w 8 -l 150 -q 5 -u 50 -n 15\
 -i "$raw"/"$sample"_1.fq.gz -I "$raw"/"$sample"_2.fq.gz\
 -o "$cleandata"/"$sample"/"$sample"_1_clean.fq.gz -O "$cleandata"/"$sample"/"$sample"_2_clean.fq.gz\
 --html "$cleandata"/"$sample"/"$sample".hmtl\
 --json "$cleandata"/"$sample"/"$sample".json > ${log}/"$sample".fastp.log 2>&1 
done


#------------------------------2_hisat2_align.sh

#!/bin/bash
for sample in `cat "$code"/sample.name`
do
"$hisat2"/hisat2 -p 10 --dta\ #--rna-strandness RF 链特异
 -x "$index"/sus_scrofa_index\
 -1 "$cleandata"/"$sample"/"$sample"_1_clean.fq.gz\
 -2 "$cleandata"/"$sample"/"$sample"_2_clean.fq.gz\
 2>"$log"/"$sample"_hisat2_align.log\
 |"$samtools"/samtools sort -@ 10\
 -O bam\
 -o "$align_output"/"$sample".align.sorted.bam - >"$log"/"$sample".samsort.log 2>&1
done

#-----------------------------3_statistic_output
#!/bin/bash
for sample in `cat ${code}/sample.name`
do
sed -n '1p' ${log}/"$sample"_hisat2_align.log|cut -d " " -f1 >> ${statistic_output}/all_clean_read.txt

sed -n '/aligned concordantly exactly 1 time/p' ${log}/"$sample"_hisat2_align.log|cut -d "(" -f1 >> ${statistic_output}/unique_reads.txt

sed -n '/aligned concordantly exactly 1 time/p' ${log}/"$sample"_hisat2_align.log|cut -d "(" -f2|cut -d ")" -f1 >> ${statistic_output}/unique_rate.txt

sed -n '/aligned concordantly >1 times/p' ${log}/"$sample"_hisat2_align.log|cut -d ")" -f1 >> ${statistic_output}/onemoretime_reads.txt
sed -n '/overall alignment rate/p' ${log}/"$sample"_hisat2_align.log|cut -d " " -f1 >> ${statistic_output}/align_rate.txt

done

paste -d "\t" \
${code}/sample.name \
${statistic_output}/all_clean_read.txt \
${statistic_output}/unique_reads.txt \
${statistic_output}/unique_rate.txt \
${statistic_output}/onemoretime_reads.txt \
${statistic_output}/align_rate.txt | sed '1i \sample.name\tall_clean_read\tunique_reads\tunique_rate\tonemoretime_reads\talign_rate' > ${statistic_output}/allsample_final_align_output_stat.txt


#-----------------------------4_stringtle_assemble.sh
#!/bin/bash
for sample in `cat "$code"/sample.name`
do
"$stringtie"/stringtie -p 6 -e\ #-rf 链特异
 -A "$assemble_output"/"$sample"_gene_abund.tab\
 -G "$ref"/Sus_scrofa.Sscrofa11.1.104.gtf\
 -o "$assemble_output"/"$sample"_stringtie.gtf\
 -l "$sample" "$align_output"/"$sample".align.sorted.bam\
 2>"$log"/"$sample"_string_assem.log;
done

#-----------------------------5_get_TPM_FTPM_rawcount.sh
#!/bin/bash
paste <(ls "$assemble_output"/*.gtf |while read line ;do basename $line _stringtie_assemble.gtf;done) <(ls "$assemble_output"/*.gtf) >"$TPM_FPKM_output"/all_sample_list.txt
python "$stringtie"/getTPM.py -i "$TPM_FPKM_output"/all_sample_list.txt -l 150 -g "$TPM_FPKM_output"/gene_tpm_matrix_nogenename.csv -t "$TPM_FPKM_output"/transcript_tpm_matrix_genename.csv 
python "$stringtie"/getFPKM.py -i "$TPM_FPKM_output"/all_sample_list.txt -l 150 -g "$TPM_FPKM_output"/gene_fpkm_matrix_nogenename.csv -t "$TPM_FPKM_output"/transcript_fpkm_matrix_genename.csv 
python ${stringtie}/prepDE.py -l 150 -i "$TPM_FPKM_output"/all_sample_list.txt -g "$TPM_FPKM_output"/gene_count_matrix.csv -t "$TPM_FPKM_output"/transcript_count_matrix.csv































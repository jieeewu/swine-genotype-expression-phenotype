library(qqman)
library(CMplot)
library(ggplot2)
library(ggforce)
library(ggprism)
library("ggrepel")
library(cowplot)
library("gridExtra")
library(dplyr)

args=commandArgs(TRUE)

tissue=args[1]
pheno_file=args[2]
predict_expr=args[3]
output_path=args[4]
test_snp_input_file=args[5]
code_path=args[6]
gene_anno_file=args[7]
eQTL_eGene_summarize=args[8]
suffix=args[9]
predict_expr_pheno <- read.table(pheno_file, header = T,check.name=F,comment.char = "")
all_predict_expr <- read.table(predict_expr, header = T)
for(pheno_line in 3:ncol(predict_expr_pheno)){
single_pheno_trio <- merge(predict_expr_pheno[,c(2,pheno_line)], all_predict_expr,by = "IID" )
pheno_name <- colnames(predict_expr_pheno)[pheno_line]
cor_result=NULL
for(gene_line in 3:ncol(single_pheno_trio)){
cor_single_pheno_trio <- cor.test(single_pheno_trio[,2],single_pheno_trio[,gene_line])
gene_name=colnames(single_pheno_trio)[gene_line]
cor_single_pheno_result <- data.frame(pheno=pheno_name,gene=gene_name,cor=cor_single_pheno_trio$estimate,pvalue=cor_single_pheno_trio$p.value, zscore = cor_single_pheno_trio$statistic)
cor_result <-rbind(cor_result,cor_single_pheno_result)
}
cor_result$FDR <- p.adjust(cor_result$pvalue, )
predict_expr_cor_result <- cor_result
write.table(cor_result, file = paste0(output_path,"/",suffix,"_",tissue,"_",pheno_name,"_association_result.txt"),quote = F, row.names = F, col.names= T, sep = "\t")

test_genotype <- read.table(test_snp_input_file,header=T,check.names = F)
test_genotype <- test_genotype[,-c(1,3:6)]
rownames(test_genotype) <- test_genotype[,1]
colnames(test_genotype) <- substr(colnames(test_genotype), 1, nchar(colnames(test_genotype)) - 2)

gene_anno <- read.table(gene_anno_file,header = T)
rownames(gene_anno) <-  gene_anno$gene_id
association_result=predict_expr_cor_result
predict_gene_expr=all_predict_expr
pheno <- predict_expr_pheno
genotype <- test_genotype
colnames(genotype) <- gsub("I","IID",colnames(genotype))


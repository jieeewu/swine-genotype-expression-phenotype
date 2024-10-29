#--------------------------gtex_tiss_chrom_training.R
argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
code_path <- argv[2]
tissue <- argv[3]
snp_input <- argv[4]
gene_input <- argv[5]
"%&%" <- function(a,b) paste(a,b, sep='')

source(code_path %&% "/15_5_gtex_v7_nested_cv_elnet.R")

snp_annot_file <- snp_input %&% "/" %&% tissue %&% "_snp_annot.chr" %&% chrom %&% ".txt"
genotype_file <- snp_input %&% "/" %&% tissue %&% "_genotype.chr" %&% chrom %&% ".txt"

gene_annot_file <- gene_input %&% "/gene_annot.parsed.txt"
expression_file <- gene_input %&% "/" %&% tissue %&% "_transformed_expression.txt"
covariates_file <- gene_input %&% "/" %&% tissue %&% "_covariates.txt"
prefix <- tissue %&% "_nested_cv"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)


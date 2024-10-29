
######################################################################
###### 
###### Calculate posterior probility
###### complete，partial，co-local and reactive
######
######################################################################

renormalize_effect_size_ratio <- function(effect_size, min_noise = 0.001) {
  ratio <- effect_size / (1 - effect_size) / (1 + min_noise)
  ratio
}

######################################################################
###### 
###### Mediation analysis
###### 
######################################################################
args <- commandArgs(TRUE)
start = as.numeric(args[1])
end = as.numeric(args[2])
tissue = args[3]
eQTL_egene_pair_file = args[4]
all_predict_expr_file = args[5]
genotype_012 = args[6]
pheno_file = args[7]
output1 = args[8]
output2 = args[9]

# Load data
eQTL_egene_pair <- read.table(eQTL_egene_pair_file, header = TRUE, check.names = FALSE)
gene_exp <- read.table(all_predict_expr_file, header = TRUE)
SNP_genotype_recode012 <- read.table(genotype_012, header = TRUE, check.names = FALSE)[, -c(1, 3:6)]
pheno <- read.table(pheno_file, header = TRUE, check.names = FALSE, comment.char = "")[, -1]

# Set rownames
rownames(gene_exp) <- gene_exp$IID
rownames(SNP_genotype_recode012) <- SNP_genotype_recode012[, 1]

library(bmediatR)
library(dplyr)
library(ggplot2)

fit_list <- list()
trio_imformation_list <- list()
mediation_result_list <- list()
mediation_post_effect_size_df <- NULL

source("calc_trio_effect_sizes.R")

for (line in start:end) {
  # Main analysis logic here...
  # Code omitted for brevity; place original analysis code here
}

save(mediation_result_list, mediation_post_effect_size_df, fit_list, trio_imformation_list, file = output1)
write.table(mediation_post_effect_size_df, file = output2, row.names = TRUE, quote = FALSE, col.names = TRUE)

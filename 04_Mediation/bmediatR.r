######################################################################
###### 
###### Calculate the proportion of variance of y explained by X
######
######################################################################
calc_effect_size <- function(y, X,
                             Z = NULL,
                             w = NULL,
                             make_intercept = TRUE) {
  
  if (is.matrix(y)) {
    y <- y[,1]
  }
  
  if (make_intercept) {
    int_mat <- matrix(rep(1, length(y)), ncol = 1); rownames(int_mat) <- names(y)
  } else{
    int_mat <- NULL
  }
  
  fit_alt <- lm.wfit(x = cbind(int_mat, Z, X), y = y, w = w)
  fit_null <- lm.wfit(x = cbind(int_mat, Z), y = y, w = w)
  
  #effect_size <- 1 - sum((y[names(fit_alt$fitted.values)] - fit_alt$fitted.values)^2)/sum((y[names(fit_alt$fitted.values)] - fit_null$fitted.values)^2)
  effect_size <- 1 - sum((y - fit_alt$fitted.values)^2)/sum((y - fit_null$fitted.values)^2)
  
  effect_size
}



######################################################################
###### 
###### Calculate empirical hyper priors for 
###### X -> M, X -> Y, and M -> Y
######
######################################################################
calc_trio_effect_sizes <- function(y, m, X,
                                   Z = NULL, Z_y = NULL, Z_m = NULL,
                                   w = NULL, w_y = NULL, w_m = NULL,
                                   make_intercept = TRUE,
                                   align_data = TRUE) {
  
  processed_data <- bmediatR:::process_data(y = y, M = m, X = X,
                                            Z_y = Z_y, Z_M = Z_m,
                                            w_y = w_y, w_M = w_m, 
                                            align_data = align_data,
                                            verbose = FALSE)
  y <- processed_data$y
  m <- processed_data$M[,1]
  X <- processed_data$X
  Z_y <- processed_data$Z_y
  Z_m <- processed_data$Z_M
  w_y <- processed_data$w_y
  w_m <- processed_data$w_M
  
  x_m <- calc_effect_size(Z = Z_m, 
                          y = m, 
                          X = X, 
                          w = w_m, 
                          make_intercept = make_intercept)
  m_y <- calc_effect_size(Z = Z_y, 
                          y = y, 
                          X = m, 
                          w = w_y, 
                          make_intercept = make_intercept)
  x_y <- calc_effect_size(Z = Z_y,
                          y = y, 
                          X = X, 
                          w = w_y, 
                          make_intercept = make_intercept)
  
  results <- c(x_m, m_y, x_y)
  names(results) <- c("x_m", "m_y", "x_y")
  results
}

renormalize_effect_size_ratio <- function (effect_size,
                                           min_noise = 0.001) {
  
  effect_size <- effect_size/(1 + min_noise)
  ratio <- effect_size/(1 - effect_size)
  ratio
}

######################################################################
###### 
###### Calculate posterior probility
###### complete，partial，co-local and reactive
######
######################################################################

args <- commandArgs(TRUE) 

genotype_012=args[1]
new_snpid=args[2]
cis_trans_both_file=args[3]
gene_exp_file=args[4]
output1=args[5]
output2=args[6]


cis_trans_both <- read.table(cis_trans_both_file)
colnames(cis_trans_both) <- c("eQTL","eGene","gene_chr","gene_TSS","SNP_chr","SNP_pos","dis","type")
cis_trans_both$dis <- abs(cis_trans_both$dis)
snp_name <- unique(cis_trans_both$eQTL)  

output_sub_rbind <- NULL

for (snp in snp_name) {
  cis_trans_both_sub <- cis_trans_both[cis_trans_both$eQTL %in% snp,]
  cis_num <- nrow(cis_trans_both_sub[which(cis_trans_both_sub[,8]=="cis"),])
  for (num in 1:cis_num){
  for (i in (cis_num + 1):nrow(cis_trans_both_sub)) {
    output_sub <- cbind(cis_trans_both_sub[num,],cis_trans_both_sub[i,])
    output_sub_rbind <- rbind(output_sub_rbind,output_sub)    
   } 
  } 
 }

output_cis_trans_pair <- output_sub_rbind[,c(1:4,7,10:12,15,16)]
colnames(output_cis_trans_pair)[c(2,6)] <- c("cis_eGene","trans_eGene")   
  
gene_exp <- read.table(gene_exp_file, comment.char = "/",sep = " ",header = T)
rownames(gene_exp) <- gene_exp$IID
SNP_genotype_recode012 <- read.table(genotype_012, header = T, check.name = F)

recode.SNP.ID <- read.table(new_snpid, header = F, check.name = F)
colnames(SNP_genotype_recode012) <- recode.SNP.ID$V1
rownames(SNP_genotype_recode012) <- SNP_genotype_recode012$FID
library(bmediatR)
library(dplyr)
trio_imformation_list <- list()
mediation_result_list <- list()
fit_list<-list()
mediation_post_effect_size_df <- NULL

rownames(output_cis_trans_pair) <- NULL
start=1
end=nrow(output_cis_trans_pair)
 for (line in start:end) {
  gap=as.numeric(line) - as.numeric(start) + 1
  eGene_exp_sub_cis <- gene_exp[ ,colnames(gene_exp) %in% output_cis_trans_pair[line,2],drop = FALSE]
  
  eGene_exp_sub_trans <- gene_exp[ ,colnames(gene_exp) %in% output_cis_trans_pair[line,6],drop = FALSE]
  eGene_exp_sub <- cbind(eGene_exp_sub_cis,eGene_exp_sub_trans)
  eGene_exp_sub$sample <- rownames(eGene_exp_sub)

  eQTL_recode012_sub <- SNP_genotype_recode012[,colnames(SNP_genotype_recode012) %in% output_cis_trans_pair[line,1], drop = F]
  eQTL_recode012_sub <- na.omit(eQTL_recode012_sub)
  eQTL_recode012_sub$sample <- rownames(eQTL_recode012_sub)
  eQTL_eGene_sub <- merge(eQTL_recode012_sub, eGene_exp_sub, by = "sample")
  eQTL_eGene_sub <- eQTL_eGene_sub[,-1]

  eQTL_name=colnames(eQTL_eGene_sub)[1]
  cis_eGene_name=colnames(eQTL_eGene_sub)[2]
  trans_eGene_name=colnames(eQTL_eGene_sub)[3]
  colnames(eQTL_eGene_sub) <- c("eQTL","cis_eGene","trans_eGene")
  trio_imformation_list[[gap]] <- eQTL_eGene_sub  
  main=paste(eQTL_name,cis_eGene_name,trans_eGene_name,sep = "_")
  names(trio_imformation_list)[gap] <- main
  trans_eGene <- subset.data.frame(eQTL_eGene_sub,select = "trans_eGene")
  cis_eGene <- subset(eQTL_eGene_sub,select = "cis_eGene")
  eQTL_012 <- subset.data.frame(eQTL_eGene_sub,select = "eQTL")
  eQTL_012_full_X <- model.matrix( ~ -1 + factor(eQTL_012$eQTL))
  true_med <- bmediatR(y = trans_eGene$trans_eGene,
                     M = cis_eGene,
                     X = eQTL_012_full_X, 
                     align_data = F, ln_prior_c = "reactive")
		
  mediation_result_list[[gap]] <- true_med  
  names(mediation_result_list)[gap] <- main

  mediation_post_df <- data.frame(eQTL=eQTL_name, cis_eGene=cis_eGene_name, trans_eGene=trans_eGene_name, complete_post=exp(true_med$ln_post_c[1,4]),
  coloc_post=exp(true_med$ln_post_c[1,7]),
  partial_post=exp(true_med$ln_post_c[1,8]),
  react_compl_post = exp(true_med$ln_post_c[1,11]),
  react_parti_post = exp(true_med$ln_post_c[1,12]))  
  rownames(mediation_post_df) <- main
  
trio_effect_sizes <- calc_trio_effect_sizes(y = trans_eGene$trans_eGene,
                                   m = cis_eGene,
                                   X = eQTL_012_full_X, 
                                   align_data = F)

trio_effect_sizes <- as.data.frame(t(trio_effect_sizes))
eQTL_eGene_sub$eQTL <- as.factor(eQTL_eGene_sub$eQTL)
correlation_Single <- cor.test(eQTL_eGene_sub$cis_eGene,eQTL_eGene_sub$trans_eGene)
correlation_Single2 <- data.frame(cor_R=correlation_Single$estimate,cor_pval=correlation_Single$p.value)
fit <- lm(eQTL_eGene_sub$cis_eGene~eQTL_eGene_sub$eQTL)
fit_list[[line]] <- fit
names(fit_list)[line] <- main

single_post_effect_size_df <- cbind.data.frame(mediation_post_df, trio_effect_sizes[1,], correlation_Single2) 
mediation_post_effect_size_df <- rbind(mediation_post_effect_size_df, single_post_effect_size_df)

plot_posterior_bar(true_med, mediator_id = "cis_eGene",
                  relabel_x = cis_eGene, main = main,
                  add_number_labels = T, label_size = 4, num_dig = 2)
  }
write.table(mediation_post_effect_size_df,file=output2,row.names=T,quote=F,col.names=T)

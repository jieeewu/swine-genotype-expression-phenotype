#---------------------------------------------MR_analysis.R
library(MendelianRandomization)
library(TwoSampleMR)
library(rlang)
library(ieugwasr)
library(MRInstruments)
library(dplyr)
library("LDlinkR")
library(Hmisc)
library("MRPRESSO")
library(stringr)

args=commandArgs(TRUE)
pheno=args[1]
tissue=args[2]
MR=args[3]
MR_output=args[4]
breed=args[5]
#-----------------------------------------------------------------------------#
#'Creates table of IVW results for the exposure-outcome combination
load(paste0(MR,"/MR_prepare/",breed,"/",breed,"_",pheno,"_",tissue,".RData"))
gene_anno_file="Sus_scrofa_11.1.104_gene_anno.txt"
genename <- read.table(paste0(tissue,"_genenum.txt"))
#gene_anno
gene_anno <- read.table(gene_anno_file,header = T)
rownames(gene_anno) <-  gene_anno$gene_id

df_result=NULL
res_sin_result=NULL
for (gene_line in 1:nrow(genename)) {
  dat=as.data.frame(MR_prepare_data[[gene_line]])
  gene_name=unlist(strsplit(names(MR_prepare_data)[gene_line],"_"))[1]
  pheno_name=unlist(strsplit(names(MR_prepare_data)[gene_line],"_"))[2]
  cat(pheno_name,tissue,gene_line,gene_name,"\n")

  dat <- dat[dat$mr_keep == TRUE, ]
  dat$beta.exposure=as.numeric(dat$beta.exposure)
  dat$se.exposure=as.numeric(dat$se.exposure)
  dat$beta.outcome=as.numeric(dat$beta.outcome)
  dat$se.outcome=as.numeric(dat$se.outcome)
  dat$eaf.exposure=as.numeric(dat$eaf.exposure)
  dat$outcome <- dat$Phenotype.outcome
  dat$exposure <- dat$gene.exposure
  dat$mr_keep <- as.logical(dat$mr_keep )

  mr_object = mr_input(bx = as.numeric(dat$beta.exposure), 
                     bxse = as.numeric(dat$se.exposure), 
                     by = as.numeric(dat$beta.outcome), 
                     byse = as.numeric(dat$se.outcome), 
                     exposure = dat$gene.exposure[1], 
                     outcome = dat$Phenotype.outcome[1], 
                     snps = dat$SNP)
#---------------------------------------------------------Sensitivity analyses

 if(nrow(dat)>1){
#Heterogeneity test
  heterogeneity_out=mr_heterogeneity(dat)
  Q_pval=heterogeneity_out[heterogeneity_out$method=="Inverse variance weighted","Q_pval"]
  if (Q_pval>0.05) {
    ivw_res_raw <- mr(dat,method_list=c('mr_ivw_fe'))
	ivw_res <- generate_odds_ratios(ivw_res_raw)
  }
  if (Q_pval<0.05) {
    ivw_res_raw <- mr(dat,method_list=c('mr_ivw_mre'))
	ivw_res <- generate_odds_ratios(ivw_res_raw)
  }
 
#leaveoneout
  mr_leaveoneout=mr_leaveoneout(dat)
  
if(is.null(ivw_res)){
  df=NA
} else {


if (ivw_res$pval < 0.001) {
  df = data.frame(exposure = ivw_res$exposure,
                  gene_symbol=gene_anno[ivw_res$exposure,"symbol"],
                  outcome = ivw_res$outcome, 
                  N_SNPs = ivw_res$nsnp,
                  Estimate = paste0(sprintf("%.2f", round(ivw_res$b, 2)), " (", sprintf("%.2f", round(ivw_res$lo_ci, 2)), ", ", sprintf("%.2f", round(ivw_res$up_ci, 2)), ")"),
                  Beta=sprintf("%.2f", round(ivw_res$b, 2)),
                  CILower=sprintf("%.2f", round(ivw_res$lo_ci, 2)),
                  CIUpper= sprintf("%.2f", round(ivw_res$up_ci, 2)),
                  P = format(ivw_res$pval, digits = 3, scientific = TRUE),
				  Q_pval = Q_pval,
				  pleiotropy_pval = mr_pleiotropy_test(dat)$pval,
				  leaveoneout_beta= mr_leaveoneout[mr_leaveoneout$SNP=="All","b"]
)
} else {
  df = data.frame(exposure = ivw_res$exposure,
                  gene_symbol=gene_anno[ivw_res$exposure,"symbol"],
                  outcome = ivw_res$outcome, 
                  N_SNPs = ivw_res$nsnp,
                  Estimate = paste0(sprintf("%.2f", round(ivw_res$b, 2)), " (", sprintf("%.2f", round(ivw_res$lo_ci, 2)), ", ", sprintf("%.2f", round(ivw_res$up_ci, 2)), ")"),
                  Beta=sprintf("%.2f", round(ivw_res$b, 2)),
                  CILower=sprintf("%.2f", round(ivw_res$lo_ci, 2)),
                  CIUpper= sprintf("%.2f", round(ivw_res$up_ci, 2)),
                  P = sprintf("%.3f", round(ivw_res$pval, 3)),
				  Q_pval = Q_pval,
				  pleiotropy_pval= mr_pleiotropy_test(dat)$pval,
				  leaveoneout_beta= mr_leaveoneout[mr_leaveoneout$SNP=="All","b"])
				  }
df_result <- rbind(df_result,df)
#----------------------------forest plot
dat$mr_keep <- as.logical(dat$mr_keep)
res_sin = mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
res_sin$exposure = gene_anno[dat$gene.exposure[1],"symbol"]
res_sin$outcome = dat$Phenotype.outcome[1]

res_sin$SNP[res_sin$SNP == "All - Inverse variance weighted (multiplicative random effects)"] = "MR-IVW estimate"
res_sin$UCL = res_sin$b + qnorm(0.975) * res_sin$se
res_sin$LCL = res_sin$b - qnorm(0.975) * res_sin$se
res_sin = res_sin[ , c("exposure", "outcome", "SNP", "b", "UCL", "LCL")]

SNPs = res_sin$SNP[!(res_sin$SNP == "MR-IVW estimate")]
SNPs_ordered = SNPs[order(res_sin$b)]
res_sin$SNP = ordered(res_sin$SNP, levels = c("MR-IVW estimate", SNPs_ordered))

res_sin_result <- rbind(res_sin_result,res_sin)

 }
} else{
ivw_res <- MendelianRandomization::mr_ivw(mr_object)
if(is.null(ivw_res)){
  df=NA
} else {
  
if (ivw_res@Pvalue < 0.001) {
  df = data.frame(exposure = ivw_res@Exposure,
                  gene_symbol=gene_anno[ivw_res@Exposure,"symbol"],
                  outcome = ivw_res@Outcome, 
                  N_SNPs = ivw_res@SNPs,
                  Estimate = paste0(sprintf("%.2f", round(ivw_res@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res@CIUpper, 2)), ")"),
                  Beta=sprintf("%.2f", round(ivw_res@Estimate, 2)),
                  CILower=sprintf("%.2f", round(ivw_res@CILower, 2)),
                  CIUpper= sprintf("%.2f", round(ivw_res@CIUpper, 2)),
                  P = format(ivw_res@Pvalue, digits = 3, scientific = TRUE),
                  Q_pval = NA,
                  pleiotropy_pval=NA,
                  leaveoneout_beta=NA)
} else {
  df = data.frame(exposure = ivw_res@Exposure, 
                  gene_symbol=gene_anno[ivw_res@Exposure,"symbol"],
                  outcome = ivw_res@Outcome, 
                  N_SNPs = ivw_res@SNPs,
                  Estimate = paste0(sprintf("%.2f", round(ivw_res@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res@CIUpper, 2)), ")"),
                  Beta=sprintf("%.2f", round(ivw_res@Estimate, 2)),
                  CILower=sprintf("%.2f", round(ivw_res@CILower, 2)),
                  CIUpper= sprintf("%.2f", round(ivw_res@CIUpper, 2)),
                  P = sprintf("%.3f", round(ivw_res@Pvalue, 3)),
                  Q_pval = NA,
                  pleiotropy_pval=NA,
                  leaveoneout_beta=NA)
				  }
df_result <- rbind(df_result,df)

#----------------------------forest plot
dat$mr_keep <- as.logical(dat$mr_keep)
res_sin = mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
res_sin$exposure = gene_anno[dat$gene.exposure[1],"symbol"]
res_sin$outcome = dat$Phenotype.outcome[1]

res_sin$SNP[res_sin$SNP == "All - Inverse variance weighted (multiplicative random effects)"] = "MR-IVW estimate"
res_sin$UCL = res_sin$b + qnorm(0.975) * res_sin$se
res_sin$LCL = res_sin$b - qnorm(0.975) * res_sin$se
res_sin = res_sin[ , c("exposure", "outcome", "SNP", "b", "UCL", "LCL")]

SNPs = res_sin$SNP[!(res_sin$SNP == "MR-IVW estimate")]
SNPs_ordered = SNPs[order(res_sin$b)]
res_sin$SNP = ordered(res_sin$SNP, levels = c("MR-IVW estimate", SNPs_ordered))

res_sin_result <- rbind(res_sin_result,res_sin)
				 }
			   }

}

save(res_sin_result,df_result,file = MR_output)

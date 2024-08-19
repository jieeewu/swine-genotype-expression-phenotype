library("data.table")
library("dplyr")
library("tidyr")
library("foreign")
library("tibble")
library("metafor")
library("meta")
library("survival")
library("ggplot2")
library("plyr")
library("gridExtra")
library("gtable")
library("grid")
library("tidyverse")
library("stringr")
library("coloc")
library("devtools")
library("glmnet")
library("MendelianRandomization")
library("TwoSampleMR")
args <- commandArgs(TRUE)
SNP_list <- read.table(args[1])
colnames(SNP_list)<-"SNP"
Outcome_exp <- read_outcome_data(
  filename = args[3], 
  snps = SNP_list$SNP, 
  sep = "\t",snp_col = "SNP",
  beta_col = "BETA",se_col = "SE",eaf_col = "AF1",
  effect_allele_col = "A1",other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",log_pval = FALSE, min_pval = 1e-200,
  chr_col = "CHR",pos_col = "POS")
  exp_dat <- read_exposure_data(
    filename = args[4],
    clump = F, sep = " ", snp_col = "SNP",
    beta_col = "BETA", se_col = "SE", 
    eaf_col = "MAF",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", 
    samplesize_col = "NMISS", min_pval = 1e-200, log_pval = FALSE, 
    chr_col = "CHR", pos_col = "BP" ,gene_col = "gene")

exp_dat_clumped <- exp_dat
exp_dat_clumped$SNP <- paste0(exp_dat_clumped$SNP,"_",tolower(exp_dat_clumped$effect_allele.exposure))
missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% Outcome_exp$SNP)]
if(length(missing_IVs) == 0) {
  print("All exposure IVs found in outcome GWAS.")
} else {
  print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
  print("List of IVs missing from outcome GWAS:")
  for (i in 1:length(missing_IVs)) {
    print(paste0(missing_IVs[i]))
  }
  missing_IVs <- toupper(missing_IVs)
  write.table(missing_IVs,file = args[5],quote = F,col.names = F,row.names = F,sep = "\n")
 }
dat <- harmonise_data(
  exposure_dat = exp_dat_clumped, 
  outcome_dat = Outcome_exp, 
  action = 2)

ld_proxies <- read.table(args[6],header = T)
fit<-try(
  Outcome_proxies <- read_outcome_data(
    filename = GWAS_file, 
    snps = ld_proxies$SNP_B, sep = "\t",snp_col = "SNP",
    beta_col = "BETA",se_col = "SE",eaf_col = "AF1",
    effect_allele_col = "A1",other_allele_col = "A2",
    pval_col = "P",
    ncase_col = "N_case",ncontrol_col = "N_control",
    samplesize_col = "N",log_pval = FALSE, min_pval = 1e-200,
    chr_col = "CHR",pos_col = "POS")
)

if("try-error" %in% class(fit)){
  print(paste0("No proxy SNP available for all exposure SNPs in outcome GWAS."))
  }else{
  print(paste0("Proxy SNP found. ", ld_proxies[ld_proxies$SNP_B %in% Outcome_proxies$SNP,"SNP_A"], " replaced with ", ld_proxies[ld_proxies$SNP_B %in% Outcome_proxies$SNP,"SNP_B"]))
  Outcome_exp <- rbind(Outcome_exp, Outcome_proxies)
  }

dat <- harmonise_data(
  exposure_dat = exp_dat_clumped, 
  outcome_dat = Outcome_exp, 
  action = 2)
  pheno_name=args[2]  
  dat <- dat[dat$mr_keep == TRUE, ]
  dat$beta.exposure=as.numeric(dat$beta.exposure)
  dat$se.exposure=as.numeric(dat$se.exposure)
  dat$beta.outcome=as.numeric(dat$beta.outcome)
  dat$se.outcome=as.numeric(dat$se.outcome)
  mr_object = mr_input(bx = as.numeric(dat$beta.exposure), 
                     bxse = as.numeric(dat$se.exposure), 
                     by = as.numeric(dat$beta.outcome), 
                     byse = as.numeric(dat$se.outcome), 
                     exposure = dat$gene.exposure[1], 
                     outcome = dat$Phenotype.outcome[1], 
                     snps = dat$SNP)
  ivw_res = MendelianRandomization::mr_ivw(mr_object) 

if (ivw_res@Outcome == pheno_name) {
  ivw_res@Estimate = exp(ivw_res@Estimate) 
  ivw_res@CILower = exp(ivw_res@CILower)
  ivw_res@CIUpper = exp(ivw_res@CIUpper)
}
if (ivw_res@Exposure == pheno_name) {
  ivw_res@Estimate = log(2)*(ivw_res@Estimate)
  ivw_res@CILower = log(2)*(ivw_res@CILower)
  ivw_res@CIUpper = log(2)*(ivw_res@CIUpper)
}
  df = data.frame(exposure = ivw_res@Exposure,
                  gene_symbol=gene_anno[ivw_res@Exposure,"symbol"],
                  outcome = ivw_res@Outcome, 
                  N_SNPs = ivw_res@SNPs,
                  Estimate = paste0(sprintf("%.2f", round(ivw_res@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res@CIUpper, 2)), ")"),
                  OR=sprintf("%.2f", round(ivw_res@Estimate, 2)),
                  CILower=sprintf("%.2f", round(ivw_res@CILower, 2)),
                  CIUpper= sprintf("%.2f", round(ivw_res@CIUpper, 2)),
                  P = format(ivw_res@Pvalue, digits = 3, scientific = TRUE))

dat$mr_keep <- as.logical(dat$mr_keep)
res_sin = mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
res_sin$exposure = gene_anno[dat$gene.exposure[1],"symbol"]
res_sin$outcome = dat$Phenotype.outcome[1]

res_sin$SNP[res_sin$SNP == "All - Inverse variance weighted (multiplicative random effects)"] = "MR-IVW estimate"
res_sin$UCL = res_sin$b + qnorm(0.975) * res_sin$se
res_sin$LCL = res_sin$b - qnorm(0.975) * res_sin$se

SNPs = res_sin$SNP[!(res_sin$SNP == "MR-IVW estimate")]
SNPs_ordered = SNPs[order(res_sin$b)]
res_sin$SNP = ordered(res_sin$SNP, levels = c("MR-IVW estimate", SNPs_ordered))

save(res_sin,df,file = args[7])
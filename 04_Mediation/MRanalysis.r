############################################################
# Mendelian randomization analysis using eQTL instruments
# IVW method
############################################################
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)
args <- commandArgs(TRUE)
SNP_list <- read.table(args[1])
colnames(SNP_list) <- "SNP"
gene_anno <- read.table(args[7], header=T, row.names=1)
trait_type <- args[8]
#---------------------- read outcome GWAS ----------------------
Outcome_exp <- read_outcome_data(
  filename=args[3],
  snps=SNP_list$SNP,
  sep="\t",
  snp_col="SNP",
  beta_col="BETA",
  se_col="SE",
  eaf_col="AF1",
  effect_allele_col="A1",
  other_allele_col="A2",
  pval_col="P",
  samplesize_col="N",
  chr_col="CHR",
  pos_col="POS",
  log_pval=FALSE,
  min_pval=1e-200)

#---------------------- read exposure eQTL ----------------------

exp_dat <- read_exposure_data(
  filename=args[4],
  clump=FALSE,
  sep=" ",
  snp_col="SNP",
  beta_col="BETA",
  se_col="SE",
  eaf_col="MAF",
  effect_allele_col="A1",
  other_allele_col="A2",
  pval_col="P",
  samplesize_col="NMISS",
  chr_col="CHR",
  pos_col="BP",
  gene_col="gene",
  log_pval=FALSE,
  min_pval=1e-200
)

exp_dat_clumped <- exp_dat
exp_dat_clumped$SNP <- paste0(exp_dat_clumped$SNP,"_", tolower(exp_dat_clumped$effect_allele.exposure))
#---------------------- check missing SNP ----------------------

missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% Outcome_exp$SNP)]
if(length(missing_IVs)>0){
  write.table(
    toupper(missing_IVs),
    file=args[5],
    quote=FALSE,
    col.names=FALSE,
    row.names=FALSE,
    sep="\n")
}
#---------------------- add proxy SNP ----------------------

if(file.exists(args[6])){
  ld_proxies <- read.table(args[6],header=T)
  if(nrow(ld_proxies)>0){
    fit <- try( Outcome_proxies <- read_outcome_data(
        filename=args[3],
        snps=ld_proxies$SNP_B,
        sep="\t",
        snp_col="SNP",
        beta_col="BETA",
        se_col="SE",
        eaf_col="AF1",
        effect_allele_col="A1",
        other_allele_col="A2",
        pval_col="P",
        samplesize_col="N",
        chr_col="CHR",
        pos_col="POS",
        log_pval=FALSE,
        min_pval=1e-200),
        silent=TRUE)
    if(!inherits(fit,"try-error")){
      Outcome_exp <- rbind(
        Outcome_exp,
        Outcome_proxies)}}}
#---------------------- harmonise ----------------------
dat <- harmonise_data(exposure_dat=exp_dat_clumped, outcome_dat=Outcome_exp, action=2)
dat <- dat[dat$mr_keep==TRUE,]
dat$beta.exposure <- as.numeric(dat$beta.exposure)
dat$se.exposure <- as.numeric(dat$se.exposure)
dat$beta.outcome <- as.numeric(dat$beta.outcome)
dat$se.outcome <- as.numeric(dat$se.outcome)
#---------------------- MR analysis ----------------------
if(length(dat$SNP)>1){
  mr_object <- mr_input(bx=dat$beta.exposure, bxse=dat$se.exposure, by=dat$beta.outcome,
    byse=dat$se.outcome,
    exposure=dat$gene.exposure[1],
    outcome=dat$Phenotype.outcome[1],
    snps=dat$SNP)

  ivw_res <- MendelianRandomization::mr_ivw(mr_object)
  Estimate <- ivw_res@Estimate
  CILower <- ivw_res@CILower
  CIUpper <- ivw_res@CIUpper
  P <- ivw_res@Pvalue
}else{
  Estimate <- dat$beta.outcome/dat$beta.exposure
  ratio_se <- sqrt(
    (dat$se.outcome/dat$beta.exposure)^2+(dat$beta.outcome*dat$se.exposure/dat$beta.exposure^2)^2)
  CILower <- Estimate-qnorm(0.975)*ratio_se
  CIUpper <- Estimate+qnorm(0.975)*ratio_se
  P <- 2*pnorm(-abs(Estimate/ratio_se))
}
if(trait_type=="binary"){
  Estimate <- exp(Estimate)
  CILower <- exp(CILower)
  CIUpper <- exp(CIUpper)}
#---------------------- result table ----------------------
gene_id <- dat$gene.exposure[1]
gene_symbol <- ifelse(gene_id %in% rownames(gene_anno), gene_anno[gene_id,"symbol"], gene_id)
df <- data.frame(
  exposure=gene_id,
  gene_symbol=gene_symbol,
  outcome=dat$Phenotype.outcome[1],
  N_SNPs=length(dat$SNP),
  Estimate=paste0(sprintf("%.2f",Estimate), " (",sprintf("%.2f",CILower), ", ",sprintf("%.2f",CIUpper),")"),OR=sprintf("%.2f",Estimate),
  CILower=sprintf("%.2f",CILower),
  CIUpper=sprintf("%.2f",CIUpper),
  P=format(P,digits=3,scientific=TRUE))

#---------------------- single SNP result ----------------------
res_sin <- mr_singlesnp(dat,all_method=c("mr_ivw_mre"))
res_sin$exposure <- gene_symbol
res_sin$outcome <- dat$Phenotype.outcome[1]
res_sin$SNP[res_sin$SNP=="All - Inverse variance weighted (multiplicative random effects)"] <- "MR-IVW estimate"
res_sin$UCL <- res_sin$b+qnorm(0.975)*res_sin$se
res_sin$LCL <- res_sin$b-qnorm(0.975)*res_sin$se
SNPs <- res_sin$SNP[res_sin$SNP!="MR-IVW estimate"]
SNPs_ordered <- SNPs[order(res_sin$b)]
res_sin$SNP <- ordered(
  res_sin$SNP,
  levels=c("MR-IVW estimate",SNPs_ordered))


save(res_sin, df, dat, file=args[9])

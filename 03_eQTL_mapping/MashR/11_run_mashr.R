#!/usr/bin/env Rscript
args <- commandArgs(TRUE) 

library(dplyr)
library(data.table)
library(ashr)
library(mashr)

file_strong <- args[1]
file_random <- args[2]
out_dir     <- args[3]
subset_n    <- as.integer(args[4])  
seed        <- as.integer(args[5])  

if (is.na(subset_n)) subset_n <- 0L
if (is.na(seed)) seed <- 202511

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- paste0(normalizePath(out_dir, mustWork = FALSE), "/")

set.seed(seed)

# ---------------------------
# Read data
# ---------------------------
z_strong <- fread(file_strong)
z_random <- fread(file_random)
z_strong$pair_id <- paste0(z_strong$SNP, ",", z_strong$gene)
z_random$pair_id <- paste0(z_random$SNP, ",", z_random$gene)

rownames(z_strong) <- z_strong$pair_id
rownames(z_random) <- z_random$pair_id

z_strong$SNP <- z_strong$gene <- z_strong$pair_id <- NULL
z_random$SNP <- z_random$gene <- z_random$pair_id <- NULL
Zs <- as.matrix(z_strong)
Zr <- as.matrix(z_random)

# ---------------------------
# Optional: subset random
# ---------------------------
if (subset_n > 0 && nrow(Zr) > subset_n) {cat("[INFO] Subsetting random to", subset_n, "rows\n")
  idx <- sort(sample(seq_len(nrow(Zr)), subset_n))
  Zr <- Zr[idx, , drop = FALSE]}

# ---------------------------
# Estimate null correlation
# ---------------------------
cat("[INFO] Estimating null correlation (Vhat)\n")
data.temp <- mash_set_data(
  Bhat = Zr,
  alpha = 1,
  zero_Bhat_Shat_reset = 1e6
)

Vhat <- estimate_null_correlation_simple(data.temp)
saveRDS(Vhat, paste0(out_dir, "Vhat.RDS"))
rm(data.temp)

# ---------------------------
# mash data objects
# ---------------------------
data.random <- mash_set_data(Bhat = Zr,alpha = 1,V = Vhat,zero_Bhat_Shat_reset = 1e6)
data.strong <- mash_set_data(Bhat = Zs,alpha = 1,V = Vhat,zero_Bhat_Shat_reset = 1e6)

# ---------------------------
# Data-driven covariances 
# ---------------------------
Ued_file <- paste0(out_dir, "U.ed.RDS")

if (file.exists(Ued_file)) {
  U.ed <- readRDS(Ued_file)
} else {
  U.pca <- cov_pca(data.strong, 4)
  U.ed  <- cov_ed(data.strong, U.pca)
  saveRDS(U.ed, Ued_file)
}

# canonical covariances
U.c <- cov_canonical(data.random)

# ---------------------------
# Fit mash model (random)
# ---------------------------
cat("[INFO] Fitting mash model on random set\n")
m.r <- mash(data.random,Ulist = c(U.ed, U.c),outputlevel = 1)

saveRDS(m.r, paste0(out_dir, "m.r.RDS"))
saveRDS(get_estimated_pi(m.r), paste0(out_dir, "estimated_pi.RDS"))

# ---------------------------
# Posterior on strong
# ---------------------------
m.s <- mash(data.strong,g = get_fitted_g(m.r),fixg = TRUE)

saveRDS(m.s, paste0(out_dir, "m.s.RDS"))
saveRDS(get_lfsr(m.s), paste0(out_dir, "lfsr.RDS"))
saveRDS(get_pm(m.s), paste0(out_dir, "pm.RDS"))
saveRDS(get_psd(m.s), paste0(out_dir, "psd.RDS"))
saveRDS(get_significant_results(m.s), paste0(out_dir, "significant_results.RDS"))
saveRDS(get_pairwise_sharing(m.s), paste0(out_dir, "pairwise_sharing.RDS"))
saveRDS(get_pairwise_sharing(m.s, factor = 0), paste0(out_dir, "pairwise_sharing_factor0.RDS"))


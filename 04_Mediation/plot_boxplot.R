library("ggplot2")

out.picture = "E:\\学习\\R-project\\YF_plan\\mediation\\bmediatR\\picture\\"
tissue = "Li"
example = "6_59259410_C_ENSSSCG00000038287_ENSSSCG00000029070"
main = paste(strsplit(example, split = "_EN")[[1]][1], 
             gene_anno[strsplit(example, split = "_")[[1]][4], 7], 
             gene_anno[strsplit(example, split = "_")[[1]][5], 7], sep = "_")

# Assume eQTL_eGene_sub is already defined from previous code
fit1 = lm(cis_eGene ~ eQTL, data = eQTL_eGene_sub)

# Get mean beta value
if (length(coef(fit1)) > 2) {
  beta1 = mean(c(coef(fit1)[[2]], coef(fit1)[[3]]))
}

# Get F-statistic p-value
summary_model = summary(fit1)
fst_pvalue1 <- format(pf(summary_model$fstatistic[1], summary_model$fstatistic[2], 
                          summary_model$fstatistic[3], lower.tail = FALSE), scientific = T, digits = 2)
if (as.numeric(fst_pvalue1) < 2.2e-16) {
  fst_pvalue1 = "< 2.2e-16"
} else {
  fst_pvalue1 = paste0("= ", fst_pvalue1)
}

# Draw boxplot
P4 <- ggplot(data = eQTL_eGene_sub, aes(eQTL, cis_eGene)) +
  geom_boxplot(aes(group = eQTL, fill = eQTL, size = eQTL), notch = T) +
  geom_jitter(aes(group = eQTL, colour = eQTL), alpha = 0.5) +
  scale_fill_manual(values = c("#669BBB", "#

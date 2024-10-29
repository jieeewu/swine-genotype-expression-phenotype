 ######################################################################
###### 
###### MS summarize
######
######################################################################
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

read_model_selection <- function(mc.selec, tissue) {
  model_select <- read.table(paste0(mc.selec, tissue, "_model_select.txt"), stringsAsFactors = FALSE)
  colnames(model_select) <- c("gene", "SNP")
  model_select %>% filter(SNP != "NA,(Intercept)")
}

calculate_gene_position <- function(eQTL_gene_melt, chr) {
  eQTL_gene_melt <- eQTL_gene_melt %>%
    mutate(gene_pos = map2_dbl(gene_chr, gene_tss, ~ ifelse(.x %in% chr$V1, .y + chr[chr$V1 == .x, "add"], NA)))
  return(eQTL_gene_melt)
}

calculate_snp_position <- function(eQTL_gene_melt, chr) {
  snp_pos <- as.data.frame(do.call(rbind, strsplit(eQTL_gene_melt$snp, "_")))
  colnames(snp_pos) <- c("snp_chr", "snp_loac")
  snp_pos <- snp_pos %>%
    mutate(snp_chr = gsub("X", "19", snp_chr),
           snp_loac = as.numeric(as.character(snp_loac)))

  eQTL_gene_melt <- cbind(eQTL_gene_melt, snp_pos) %>%
    mutate(snp_pos = map2_dbl(snp_chr, snp_loac, ~ ifelse(.x %in% chr$V1, .y + chr[chr$V1 == .x, "add"], NA)))

  return(eQTL_gene_melt)
}

args = commandArgs(TRUE)
tiss = c("Mu", "Li", "AF", "BF")
Tissue2 = c("LD Muscle", "Liver", "Abdominal Adipose", "Backfat")

mc.selec = args[1]
gene.tss.tes.file = args[2]
chr.file = args[3]
out.picture = args[4]
out.file = args[5]

chr <- read.table(chr.file, header = TRUE)
chr <- chr %>%
  mutate(add = c(0, cumsum(head(V3, -1))),
         mid = V3 / 2 + add)

# Loop through tissues
summarize_output <- lapply(1:4, function(ti) {
  tissue <- tiss[ti]
  tissue2 <- Tissue2[ti]
  
  model_select <- read_model_selection(mc.selec, tissue)
  snp <- strsplit(model_select$SNP, split = ",")
  names(snp) <- model_select$gene
  
  # Melt SNP-gene relationship
  eQTL_gene_melt <- melt(snp) %>%
    rename(snp = Var1, gene = Var2) %>%
    mutate(snp = gsub("VAR_", "", snp))

  # Read TSS and TES data
  gene.tss <- read.table(gene.tss.tes.file, header = FALSE, as.is = TRUE, row.names = 1) %>%
    setNames(c("CHR", "Chain", "TSS", "TES")) %>%
    mutate(CHR = gsub("X", "19", CHR))
	
  eQTL_gene_melt <- calculate_gene_position(eQTL_gene_melt, gene.tss)
  eQTL_gene_melt <- calculate_snp_position(eQTL_gene_melt, chr)
  
  # Annotation of cis and trans
  eQTL_gene_melt <- eQTL_gene_melt %>%
    mutate(cis_trans = ifelse(gene_chr == snp_chr, 
                               ifelse(abs(snp_loac - gene_tss) <= 1e6, "cis", "trans"), 
                               "trans"),
           same_CHR = ifelse(gene_chr == snp_chr, "yes", "no")) %>%
    select(-c(gene_chr, snp_chr)) %>%
    na.omit()

  # Summarize results
  summarize <- c(tissue = tissue,
                 eGene = length(unique(eQTL_gene_melt$gene)),
                 eQTL = length(eQTL_gene_melt$snp),
                 table(eQTL_gene_melt$cis_trans),
                 same_chr = table(eQTL_gene_melt$same_CHR)[[1]])

  return(summarize)
})

summarize_output_df <- do.call(rbind, summarize_output)

# Plotting
p3 <- ggplot(eQTL_gene_melt, aes(x = snp.pos, y = gene.pos, colour = cis_trans)) +
  geom_point(size = 0.0001) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black", size = 24, family = "sans", face = "bold"),
    axis.text.y = element_text(color = "black", size = 24, family = "sans", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 45, hjust = 0.45, vjust = 0.8, family = "sans", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  ) +
  scale_color_manual(values = c("#BB0021FF", "#3C5488FF")) +
  scale_x_continuous(breaks = chr$mid, labels = c(seq(1, 18), "X")) +
  scale_y_continuous(breaks = chr$mid, labels = c(seq(1, 18), "X")) +
  labs(x = "", y = "", title = tissue2)

eQTL_gene_melt$eQTL_clust <- cut(as.numeric(eQTL_gene_melt$snp.pos), breaks = 100)
eQTL_clust <- as.data.frame(table(eQTL_gene_melt$eQTL_clust)) %>%
  mutate(label = seq(1, 100, 1))

p4 <- ggplot(eQTL_clust, aes(x = label, y = Freq)) +
  geom_line(color = "#631879FF", linewidth = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 16, family = "sans", face = "bold"),
        axis.title.x = element_text(size = 35, family = "sans", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  xlab("eQTL (chr)") +
  ylab("")

eQTL_gene_melt$gene_clust <- cut(as.numeric(eQTL_gene_melt$gene.pos), breaks = 100)
gene_clust <- as.data.frame(table(eQTL_gene_melt$gene_clust)) %>%
  mutate(label = seq(1, 100, 1))

p5 <- ggplot(gene_clust, aes(x = label, y = -Freq)) +
  geom_line(color = "#631879FF", linewidth = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 16, angle = -90, family = "sans", face = "bold"),
        axis.title.y = element_text(size = 35, family = "sans", face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  xlab("eGene-TSS (chr)") +
  ylab("") +
  scale_y_continuous(limits = c(-160, 0), breaks = seq(-150, 0, 50), labels = c(150, 100, 50, 0)) +
  scale_x_continuous(breaks = NULL) +
  coord_flip()

# Save plots 
pdf(paste0(out.picture, tissue, "_eGene_eQTL.pdf"), width = 15.5, height = 14)
print(p5 + p3 + plot_spacer() + p4 + plot_layout(widths = c(0.5, 10), heights = c(10,

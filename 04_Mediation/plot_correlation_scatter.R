# Load necessary libraries
library("bmediatR")
library("ggplot2")

# Set output directory and tissue type
out.picture <- "E:\\学习\\R-project\\YF_plan\\mediation\\bmediatR\\picture\\"
tissue <- "Li"
example <- "6_59259410_C_ENSSSCG00000038287_ENSSSCG00000029070"

# Extract gene information from example
cis_eGene <- strsplit(example, split = "_")[[1]][4]
trans_eGene <- strsplit(example, split = "_")[[1]][5]
eQTL <- strsplit(example, split = "_EN")[[1]][1]

# Load data
load(paste0(tissue, "_mediation_post_effect_size_react.RData"))

# Get eQTL information
eQTL_eGene_sub <- trio_imformation_list[[example]]
correlation_Single <- mediation_post_effect_size_df[example, ]
eQTL_eGene_sub$eQTL <- as.factor(eQTL_eGene_sub$eQTL)
pvalue <- format(correlation_Single$cor_pval, scientific = TRUE, digits = 2)

# Create scatter plot and save to PDF
pdf(paste0(out.picture, tissue, "_", main, "_all_picture_colocal2.pdf"), width = 5, height = 5)
P3 <- ggplot(aes(x = cis_eGene, y = trans_eGene, colour = eQTL), data = eQTL_eGene_sub) +
  geom_point(alpha = 0.7, shape = 20) +
  scale_color_manual(values = c("#669BBB", "#993333", "#E5C494")) +
  stat_smooth(method = lm, colour = "black") +
  xlab(cis_eGene_name) +
  ylab(trans_eGene_name) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA),
        axis.line.y.left = element_line(colour = "black", linewidth = 1.2),
        axis.line.x.bottom = element_line(colour = "black", linewidth = 1.2),
        text = element_text(size = 30),
        legend.position = "none",
        axis.text = element_text(family = "sans", face = 1, size = 20, colour = "black"),
        axis.title = element_text(family = "sans", face = 1, size = 28)
  ) +
  annotate("text", family = "sans", size = 6, label = paste("R =", round(correlation_Single$cor_R, digits = 2)),
           x = -1.8, y = max(eQTL_eGene_sub$cis_eGene, eQTL_eGene_sub$trans_eGene)) +
  annotate("text", family = "sans", size = 6, label = paste("P =", pvalue),
           x = -1.8, y = max(eQTL_eGene_sub$cis_eGene, eQTL_eGene_sub$trans_eGene) - 0.25)
print(P3)
dev.off()

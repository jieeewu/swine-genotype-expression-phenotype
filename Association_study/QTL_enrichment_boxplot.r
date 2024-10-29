library(ggplot2)
library(dplyr)
library(ggrepel)
cbbPalette <- c("#669BBB", "#6A2C70", "#F08A5D", "#5CB85CFF", "#A6761D")

input_file <- read.table("enrichment_in_QTL.txt", header = TRUE)  # Ensure to read the input again
input_file$tissue <- as.factor(input_file$tissue)
input_file$tissue <- factor(input_file$tissue, levels = c("Mu", "AF", "Li", "BF"))
input_file$group <- as.factor(input_file$group)

# Annotate
desired_FC <- 2
data <- input_file %>%
  mutate(is_annotate = ifelse(mean_ratio > desired_FC, "yes", "no"))
data_sub <- subset(data, is_annotate == "yes")
data$is_annotate <- as.factor(data$is_annotate)
data_sub <- data_sub %>%
  mutate(pheno = case_when(
    pheno == "BodyHeight_cm" ~ "BodyH",
    pheno == "CarcassStrL_cm" ~ "CarLen",
    pheno == "LD_ColorM_a12h" ~ "a*12hr",
    pheno == "CC_cm" ~ "ChestC",
    pheno == "LD_ColorScore_12h" ~ "MCS12hr",
    pheno == "LD_ColorScore_45min" ~ "MCS45min",
    pheno == "RibBackFat_cm" ~ "RibFT",
    pheno == "ThoracicLength_cm" ~ "ThoraVLen",
    pheno == "WaistBackFat_cm" ~ "WaistFT",
    TRUE ~ pheno
  ))

data$group <- factor(data$group, levels = c("Carcass", "Exterior", "MeatQuality", "Reproduction", "Production"))

pdf("E:\\学习\\R-project\\YF_plan\\GWAS\\eQTL_enrich_in_QTL\\QTL_enrichment_all_pheno_boxplot.pdf", width = 4, height = 6)
ggplot() +
  geom_boxplot(data = data, aes(tissue, mean_ratio, group = tissue), notch = FALSE) +
  geom_jitter(data = data, aes(tissue, mean_ratio, group = tissue, color = group), alpha = 0.8) +
  scale_color_manual(values = cbbPalette, labels = c("Carcass", "Exterior", "Meat quality", "Reproduction", "Production")) +
  xlab("") +
  ylab("Fold enrichment") +
  ylim(0, 2.6) +
  scale_x_discrete(labels = c("LD Muscle", "Abdominal Adipose", "Liver", "Backfat")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA),
        axis.line.y.left = element_line(colour = "black", linewidth = 1.2),
        axis.line.x.bottom = element_line(colour = "black", linewidth = 1.2),
        text = element_text(size = 15, family = "sans", face = 1),
        legend.position = c(0.5, 1.03),
        legend.title = element_blank(),
        legend.key.width = unit(0.2, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.text = element_text(family = "sans", colour = "black", face = 1, size = 15),
        axis.text.x = element_text(family = "sans", angle = 30, colour = "black", face = 1, size = 15, vjust = 0.98, hjust = 0.95),
        axis.text.y = element_text(family = "sans", colour = "black", face = 1, size = 15),
        axis.title = element_text(family = "sans", face = 1, size = 20),
        plot.margin = margin(t = 1, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) +
  guides(color = guide_legend(nrow = 2)) +
  geom_text_repel(data = data_sub, aes(x = tissue, y = mean_ratio, label = pheno, color = group),
                  size = 3.3,
                  force_pull = 0,
                  nudge_y = 0.05,
                  direction = "x",
                  show.legend = FALSE,
                  segment.size = 0.2,
                  max.iter = 10e4, max.time = 0.5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 5))

dev.off()

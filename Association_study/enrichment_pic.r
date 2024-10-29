#===============================================================================================================
#==================================== enrichment of eQTLs in QTL ===============================================
#===============================================================================================================
library(ggplot2)
reult_merge <- read.table("enrichment_in_QTL.txt", header = TRUE)
reult_merge$pheno <- as.factor(reult_merge$pheno)
cbbPalette <- c("#C11221", "#669BBB", "#6A2C70", "#F08A5D", "#5CB85CFF")
pheno_name <- read.table("pheno_name.txt", header = TRUE, sep = "\t")
reult_merge$pheno <- factor(reult_merge$pheno, levels = pheno_name$pheno1)
rownames(pheno_name) <- pheno_name$pheno1
reult_merge$group <- pheno_name[reult_merge$pheno, "group"]
reult_merge$pheno2 <- pheno_name[reult_merge$pheno, 3]
reult_merge$pheno2 <- factor(reult_merge$pheno2, levels = pheno_name$pheno3)

lable <- function(x) { stringr::str_wrap(x, width = 35) }
input_file <- reult_merge

#-------------------------------plot
pdf("QTL_enrichment.pdf", width = 15, height = 6)
p2 <- ggplot() +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(family = "sans", face = 1, colour = 'black', size = 15),
        axis.text.y = element_text(colour = 'black', size = 15, family = "sans", face = 1),
        plot.margin = margin(t = 1, r = 1, b = 3, l = 1, unit = "cm"),
        axis.text.x = element_text(colour = 'black', family = "sans", face = 1, angle = 30, size = 10, vjust = 0.98, hjust = 0.95),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key = element_blank(),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 15, family = "sans", face = 1, colour = "black", hjust = 0.5)) +
  scale_x_discrete(limits = levels(input_file$pheno), labels = lable(levels(input_file$pheno2))) +
  xlab("") +
  ylab("Fold enrichment") +
  labs(title = "")
for (i in 1:(length(unique(input_file$pheno)) - 1)) {
  p2 <- p2 + annotate('rect', xmin = i + 0.5, xmax = i + 1.5, ymin = -Inf, ymax = Inf,
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
}
p2 <- p2 +
  geom_point(data = input_file, aes(x = pheno, y = mean_ratio, fill = tissue, colour = tissue),
             position = position_dodge(0.8), shape = 21, size = 3, alpha = 1) + 
  geom_hline(aes(yintercept = 1), colour = "grey", linetype = "dashed") +
  scale_y_continuous(limits = c(0, 3.5), breaks = c(0:4)) +
  scale_fill_manual(values = cbbPalette, labels = c("Abdominal Adipose", "Backfat", "Liver", "LD Muscle", "All")) +
  scale_color_manual(values = cbbPalett

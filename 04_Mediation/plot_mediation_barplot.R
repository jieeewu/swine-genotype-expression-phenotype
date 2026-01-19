# Load necessary libraries
library("bmediatR")
library("ggplot2")
library(RColorBrewer)

tissue <- "Li"
example <- "6_59259410_C_ENSSSCG00000038287_ENSSSCG00000029070"

# Extract gene information from example
cis_eGene <- strsplit(example, split = "_")[[1]][4]
trans_eGene <- strsplit(example, split = "_")[[1]][5]
eQTL <- strsplit(example, split = "_EN")[[1]][1]

# Load data
load(paste0(tissue, "_mediation_post_effect_size_react.RData"))
gene_anno <- read.table("Sus_scrofa_11.1.104_gene_anno_updata.txt", header = TRUE)
rownames(gene_anno) <- gene_anno$gene_id

# Get mediation results
true_med <- mediation_result_list[[example]]
cis_eGene_name <- gene_anno[cis_eGene, 7]
trans_eGene_name <- gene_anno[trans_eGene, 7]
main <- paste(eQTL, cis_eGene_name, trans_eGene_name, sep = "_")

# Prepare mediation post dataframe
mediation_post_df <- data.frame(
  cis_eGene = cis_eGene_name,
  trans_eGene = trans_eGene_name,
  complete_post = exp(true_med$ln_post_c[1, 4]),
  partial_post = exp(true_med$ln_post_c[1, 8]),
  coloc_post = exp(true_med$ln_post_c[1, 7]),
  react_compl_post = exp(true_med$ln_post_c[1, 11]),
  react_parti_post = exp(true_med$ln_post_c[1, 12])
)
mediation_post_df$other <- 1 - sum(mediation_post_df[3:7])
mediation_post_df <- as.data.frame(t(mediation_post_df))
mediation_post_df[3:nrow(mediation_post_df), ] <- round(as.numeric(mediation_post_df[3:nrow(mediation_post_df), ]), digits = 2)

# Define colors
color <- c("#669BBB", "#C7B4D9", "#002F49", "#FEF0D5", "#E5C494", "#B3B3B3")
fill <- "black"
max_lim <- max(as.numeric(mediation_post_df[3:nrow(mediation_post_df), ]))

# Create bar plot and save to PDF
pdf(paste0(out.picture, tissue, "_", main, "_all_picture_colocal1.pdf"), width = 6, height = 5)
par(omi = c(0.5, 0.5, 0.5, 0.5), mar = c(4, 4, 4, 10), family = "sans", font = 1)
P2 <- barplot(height = as.numeric(mediation_post_df[3:nrow(mediation_post_df), ]),
              family = "sans", las = 1, ylim = c(0, max_lim + 0.1),
              cex.axis = 1, family = "sans", font = 1,
              col = color, border = fill)
box(lty = 1, lwd = 2, bty = "l")
legend(x = 6.8, y = max_lim,
       legend = c("complete med", "partial med", "co-local", "complete med (react)", "partial med (react)", "other non-med"),
       fill = color, border = fill,
       bty = "n", cex = 0.7, text.width = 0.2, xpd = TRUE)
text(x = P2, y = as.numeric(mediation_post_df[3:nrow(mediation_post_df), ]) + 0.05,
     labels = as.numeric(mediation_post_df[3:nrow(mediation_post_df), ]),
     family = "sans", font = 1)
axis(side = 1, las = 1, cex.axis = 1.3, family = "sans", at = 3.5, labels = cis_eGene_name,
     font = 1, line = 0, tick = FALSE)
dev.off()

print(P2)

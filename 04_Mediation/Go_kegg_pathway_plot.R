library(dplyr)
library(ggplot2)
library(stringr)
library(tools)
library(ggtext)   
library(cowplot)
library(grid)

df <- read.table("Enrichment_result.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

df$Category <- gsub("KEGG_PATHWAY", "KEGG", df$Category)
df$Category <- gsub("GOTERM_MF_DIRECT", "GO_MF", df$Category)
df$Category <- gsub("GOTERM_CC_DIRECT", "GO_CC", df$Category)
df$Category <- gsub("GOTERM_BP_DIRECT", "GO_BP", df$Category)


df <- df %>%
  mutate(
    Term = gsub("^ssc[0-9]+:", "", Term),
    Term = gsub("^GO:[0-9]+~", "", Term),
    Term = gsub("^([a-z])", "\\U\\1", Term, perl = TRUE),
	Term = stringr::str_wrap(Term, width = 48),
    Term = gsub("\n", "<br>", Term),
    
    Rich_factor = Count / `List.Total`
  )


tissues <- unique(df$Tissue)

cat_colors <- c(
  "KEGG"  = "#1f78b4",
  "GO_MF" = "#33a02c",
  "GO_CC" = "#6a3d9a",
  "GO_BP" = "#ff7f00"
)

cat_line_colors <- c(
  "KEGG"  = alpha("#1f78b4", 0.25),
  "GO_MF" = alpha("#33a02c", 0.25),
  "GO_CC" = alpha("#6a3d9a", 0.25),
  "GO_BP" = alpha("#ff7f00", 0.25)
)

cat_order <- c("KEGG", "GO_MF", "GO_CC", "GO_BP")

# ============================================================================
for (t in tissues) {
  
  df_t <- df %>%
    filter(Tissue == t) %>%
    mutate(Category = factor(Category, levels = cat_order)) %>%
    arrange(Category, PValue) %>%
    mutate(Term = factor(Term, levels = rev(unique(Term))))
  
  label_tbl <- df_t %>% distinct(Term, Category)
  term_label_fun <- function(terms) {
    sapply(terms, function(tt){
      cat <- label_tbl$Category[label_tbl$Term == tt][1]
      col <- cat_colors[cat]
      paste0("<span style='color:", col, "'>", tt, "</span>")
    })
  }
  
  df_lines <- df_t %>%
    mutate(
      y_pos = as.numeric(Term),
      line_col = cat_line_colors[Category]
    )
  
  p <- ggplot(df_t,
              aes(x = -log10(PValue),
                  y = Term,
                  size = Count,
                  color = Rich_factor)) +
    
  geom_point(alpha = 0.9) +
    
    scale_size_continuous(range = c(3, 10)) +
    scale_color_gradient(low = "#085CA0", high = "#B1182C", name = "Rich factor") +
    
    scale_y_discrete(labels = term_label_fun) +
    theme_bw() +
    
    labs(
      title = t,
      x = expression(-log[10](italic(P)~value)),
      y = "",
      size = "Gene count"
    ) +
    
    theme(
      text = element_text(family = "sans", size = 12),
      plot.title = element_text(family = "sans", size = 14, face = "bold"),
      axis.text.y = ggtext::element_markdown(size = 10),
      axis.text = element_text(family = "sans", size = 10),
      legend.text = element_text(family = "sans", size = 10),
      legend.title = element_text(family = "sans", size = 11),
      legend.position = "right",
      plot.margin = margin(10, 60, 10, 10)
    ) +
    
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.2)))
 
  outfile_pdf <- paste0("Figure_", gsub(" ", "_", t), ".pdf")
  outfile_png <- paste0("Figure_", gsub(" ", "_", t), ".png")
  
  ggsave(outfile_pdf, p, width = 7, height = 6, dpi = 300)
  ggsave(outfile_png, p, width = 7, height = 6, dpi = 300)

}















library("grid");
library(VennDiagram);
library(RColorBrewer);

pic.out="/mnt/research/qgg/wu/YFJH/eQTL_mapping/eQTL_output/FDR/common_gene/"
for(tissue in c("Muscle","Liver")){
dataset1=read.table(paste0(,tissue,"_gene_Ballester.txt"))
dataset2=read.table(paste0(,tissue,"_gene_pigGTEx.txt"))
SGEP=read.table(paste0(tissue,"_gene_SGEP.txt"))


venn_list <- list("SGEP"=SGEP$V1,"Ballester et al."=dataset1$V1, "pigGTEx"=dataset2$V1)

pdf(paste0(pic.out, tissue,"_venn_diagram_according_to_genes.pdf"),width = 8,height = 8)

P_veen<- venn.diagram(venn_list,
                      fill = c("#199FB1","#5CB85CFF","#D43F3AFF"),
                      alpha = c(0.5,0.5,0.5),
                      cex = 2.5,
                      cat.fontfamily = "sans",
                      cat.fontface=2,
                      fontfamily="sans",
                      fontface=2,
                      cat.cex = 2.5,
					  main=tissue,
					  main.cex=3,
					  main.fontfamily = "sans",
					  main.fontface=2,
                      filename = NULL)

grid.draw(P_veen)

dev.off()

}

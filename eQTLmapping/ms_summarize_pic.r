 ######################################################################
###### 
###### MS summarize
######
######################################################################
args = commandArgs(TRUE)
tiss=c("Mu","Li","AF","BF")
Tissue2=c("LD Muscle","Liver","Abdominal Adipose","Backfat")

mc.selec=args[1]
gene.tss.tes.file=args[2]
chr.file=args[3]
out.picture=args[4]
out.file=args[5]

#input
single_snp_gene=data.frame(snp.number=1:6)
all_snp_gene=NULL
for (ti in 1:4){
tissue=tiss[ti]
tissue2=Tissue2[ti]
summarize_output=NULL
model_select <- read.table(paste0(mc.selec,tissue,"_model_select.txt"))
model_select <- apply(model_select, 2, as.character)
colnames(model_select) <- c("gene","SNP")
model_select <- model_select[which(model_select[,2]!="NA,(Intercept)"), ]
snp<- strsplit(model_select[, 2], split = ",")
names(snp) <- model_select[,1]

#snp-eGene
library(reshape2)
eQTL_gene_melt <- melt(snp)
colnames(eQTL_gene_melt) <- c("snp","gene")
eQTL_gene_melt$snp <- gsub("VAR_", "", eQTL_gene_melt$snp)

#tss tes
gene.tss <- read.table(gene.tss.tes.file, header = FALSE, as.is = TRUE, row.names = 1) #gene, CHR, +/-, start_pos，end_pos
colnames(gene.tss) <- c("CHR", "Chain", "TSS", "TES")
gene.tss$CHR <- gsub("X", "19", gene.tss$CHR)

#chr
chr <- read.table(chr.file, header = T)
chr[1,4]=0
chr[1,5] <- chr[1,3]/2 + chr[1,4] 
#chr[1,6]=chr[1,5]
for (k in 2:nrow(chr)) {
  chr[k,4] <- chr[k-1,3] + chr[k-1,4] 
  chr[k,5] <- chr[k,3]/2 + chr[k,4] 
}
colnames(chr)[4:5] <- c("add","mid")
#gene TSS
eQTL_gene_melt[,c("gene_chr","gene_tss")] <- gene.tss[eQTL_gene_melt$gene,c(1,3)]


for (m in 1:19) {
  for (n in 1:nrow(eQTL_gene_melt)) {
    if(eQTL_gene_melt[n, 3]==chr[m,1]){
      eQTL_gene_melt[n, 5] = eQTL_gene_melt[n, 4] + chr[m,4] 
    }
  }
}
colnames(eQTL_gene_melt)[5] <- "gene.pos"
eQTL_gene_melt<- eQTL_gene_melt[which(eQTL_gene_melt$gene_chr %in% seq(1,19,1)),]
eQTL_gene_melt$gene_chr <- as.numeric(eQTL_gene_melt$gene_chr)
eQTL_gene_melt <- na.omit(eQTL_gene_melt)
rownames(eQTL_gene_melt)<-NULL

#snp
snp_pos <- as.data.frame(do.call(rbind, strsplit(eQTL_gene_melt$snp,"_" )))
snp_pos[,1] <- gsub("X", "19", snp_pos[,1])
snp_pos[,2] <- as.numeric(as.character(snp_pos[,2]))
eQTL_gene_melt[ ,c("snp_chr", "snp_loac")] <- snp_pos[ ,1:2]
for (m in 1:19) {
  for (n in 1:nrow(eQTL_gene_melt)) {
    if(eQTL_gene_melt[n, 6] == chr[m,1]){
      eQTL_gene_melt[n, 8] = eQTL_gene_melt[n, 7] + chr[m,4]  #+(m-1)*3000000
    }
  }
}

colnames(eQTL_gene_melt)[8] <- "snp.pos"

dim(eQTL_gene_melt)

#anno cis and trans
eQTL_gene_melt[,3:8] <- apply(eQTL_gene_melt[,3:8],2,as.numeric)
eQTL_gene_pos <- eQTL_gene_melt

#i=6
for (i in 1:nrow(eQTL_gene_melt)) {
  
  if (eQTL_gene_melt[i,3]==eQTL_gene_melt[i, 6]) {
    tss=as.numeric(eQTL_gene_melt[i, 4])
    eQTL_gene_pos[i ,9] <- NA
    eQTL_gene_pos[i ,10] <- "yes"
    if (eQTL_gene_melt[i, 7] <= tss + 1000000 & eQTL_gene_melt[i, 7] >= tss - 1000000) {
      eQTL_gene_pos[i, 9] <- "cis"
    } else{
      eQTL_gene_pos[i, 9] <- "trans"
    }
  } else{
    eQTL_gene_pos[i, 9] <- "trans"
  }
}

colnames(eQTL_gene_pos)[9:10] <- c("cis_trans", "same_CHR")
table(eQTL_gene_pos[, 9] )
eQTL_gene_pos[, 9] <- as.factor(eQTL_gene_pos[, 9])

summarize <- t(as.matrix(c(tissue=tissue,eGene=length(unique(eQTL_gene_pos$gene)),eQTL=length(eQTL_gene_pos$snp),table(eQTL_gene_pos$cis_trans),same_chr=table(eQTL_gene_pos$same_CHR)[[1]])))

library(ggplot2)
p3 <- ggplot(data = eQTL_gene_pos, aes(x = snp.pos, y = gene.pos, colour = cis_trans))+ #  shape = cis_trans
  geom_point(size=0.0001)+
  theme_bw()+
  theme(
        axis.text.x = element_text(color = "black", size = 24,family = "sans",face = "bold"),
        axis.text.y = element_text(color = "black", size = 24,family = "sans",face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
		plot.title = element_text(size = 45,hjust = 0.45, vjust = 0.8,family = "sans",face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none"  
        ) +   
  scale_color_manual(values = c("#BB0021FF","#3C5488FF"))+
  scale_x_continuous(breaks = chr[1:19,5],
                     labels = c(seq(1,18,1), "X"),
                     )+ 
  scale_y_continuous(breaks = chr[1:19,5],
                     labels = c(seq(1,18,1), "X"),
                     )+
  xlab("")+
  ylab("")+
  labs(title = tissue2)


eQTL_gene_pos$eQTL_clust <- cut(as.numeric(eQTL_gene_pos$snp.pos), breaks = 100)
eQTL_clust <- as.data.frame(table(eQTL_gene_pos$eQTL_clust))
eQTL_clust$label <- seq(1,100,1)
p4 <- ggplot(eQTL_clust, aes(x = label, y = Freq))+
  geom_line(color="#631879FF", linewidth=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
		axis.text.y = element_text(color = "black",size = 16,family = "sans",face = "bold"),
		axis.title.x = element_text(size = 35,family = "sans",face = "bold"),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
		)+
  scale_x_continuous(breaks = NULL)+
  xlab("eQTL (chr)")+
  ylab("")


eQTL_gene_pos$gene_clust <- cut(as.numeric(eQTL_gene_pos$gene.pos), breaks = 100)
gene_clust <- as.data.frame(table(eQTL_gene_pos$gene_clust))
gene_clust$label <- seq(1,100,1)

p5 <- ggplot(gene_clust, aes(x = label, y = -Freq))+
  geom_line(color="#631879FF", linewidth=0.5)+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",size = 16,angle = -90,family = "sans",face = "bold"),
		axis.title.y = element_text(size = 35,family = "sans",face = "bold"),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
		)+
  xlab("eGene-TSS (chr)")+
  ylab(" ")+
  scale_y_continuous(limits = c(-160,0),breaks = seq(-150,0,50), labels = c(150,100,50,0))+
  scale_x_continuous(breaks = NULL)+
  coord_flip()  

library("ggpubr")
library(patchwork)
pdf(paste0(out.picture,tissue,"_eGene_eQTL.pdf"), width= 15.5, height = 14)
print(p5+p3+plot_spacer()+p4+plot_layout(widths=c(0.5,10),heights = c(10,0.5),nrow = 2,ncol=2))
dev.off()















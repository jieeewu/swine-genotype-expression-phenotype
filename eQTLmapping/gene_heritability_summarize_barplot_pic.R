#-------------------------------------------------------------------------------------------------------------------------------
#---------------------------------plot heritability summarize bar plot ---------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
args = commandArgs(TRUE)
Tissue<-c("LD Muscle","Liver","Abdominal Adipose","Backfat","Mu","Li","AF","BF")
lable=c("a","b","c","d")
pdf(paste0(args[1],"/","Gene_heritability_hist.pdf"),width = 18,height = 5)
par(mfrow = c(1,4),omi=c(0.5,0.5,0.5,0.5),mar=c(4,4,4,4),family = "sans",font = 1)
for (i in 1:4){
tissue=Tissue[i]
Tiss=Tissue[i+4]
heri_qcov_sva_sex <- read.table(paste0(args[2],"/",Tiss,"_gene_heritability.txt"), header = F)
colnames(heri_qcov_sva_sex) <- c("gene", "h2", "SE", "pval")
heri_qcov_sva_sex$FDR <- p.adjust(heri_qcov_sva_sex$pval, method = "BH")
heri_qcov_sva_sex_FDR <- heri_qcov_sva_sex[which(heri_qcov_sva_sex$FDR<0.05),]
write.table(heri_qcov_sva_sex_FDR$gene,paste0(args[2],"/",Tiss,"_genename_heri_FDR05.txt"),quote = F,col.names = F,row.names = F)

hist(heri_qcov_sva_sex$h2, xlab = "", ylab = "", family = "sans",ylim = c(0,4000),
     font = 1, cex.lab = 2.5, las = 1, cex.axis = 2,cex.main=2.5,font.main=1, main = tissue, col="white", border="black")
hist(heri_qcov_sva_sex_FDR$h2, add=T, col = "#C11221") #BB0021FF
legend("topright",legend = "FDR = 0.05", fill = "#C11221", bty = "n", cex=2)
mtext("Heritability",side = 1, line = 4, adj = 0.5,cex = 1.5)
mtext("Count",side = 2, line = 5, adj = 0.5,cex = 1.5)
box(lty = 1, lwd = 2,bty="l")
#fig_label(lable[i],pos="topleft", cex=2) 
}
dev.off()

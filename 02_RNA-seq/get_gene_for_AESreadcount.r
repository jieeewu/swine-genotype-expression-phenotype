
tpm<-read.csv("gene_tpm_matrix.csv",row.names=1)
LD_tpm_mean <- apply(LD_tpm,1,mean)
tpm_mean2<-LD_tpm_mean[!is.na(LD_tpm_mean)]
tpm_mean2 <- as.data.frame(tpm_mean2)
write.table(tpm_mean2,"tpm_mean.txt",col.names=F,row.names=T)
length(tpm_mean2[tpm_mean2>=200,])

library("ggplot2")
colnames(tpm_mean2)
tpm_mean2$gene_id <- rownames(tpm_mean2)
pdf("tpm_mean_point.pdf")
ggplot()+
  geom_point(data=tpm_mean2,aes(x=gene_id,y=tpm_mean2),shape=21,size=1)+
  geom_hline(aes(yintercept=c(200,400)),colour="red",linetype=2)
dev.off()

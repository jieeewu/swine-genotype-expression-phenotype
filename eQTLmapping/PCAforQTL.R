# =========================================================
# = Gene Expression Adjust Analysis =
# =========================================================
args = commandArgs(TRUE)
tissue = args[1]
tpm <- read.csv(args[2],header = T,row.names = 1, check.names = F)
load(args[3])
group <- get(paste0(tissue,"_group"))
sample_keep <- tpm[,which(colnames(tpm) %in% as.matrix(group$Sample_id))]
unique(colnames(sample_keep)==group$Sample_id)
tpm_filter <- sample_keep[rowSums(sample_keep>1)>=0.75*ncol(sample_keep),]
genename <- rownames(tpm_filter)
write.table(genename, file = args[4], quote = F, col.names = F,row.names = F)

tpm_filter_mean <- as.data.frame(apply(tpm_filter, 1, mean))
colnames(tpm_filter_mean) <- c("mean_tpm")

pdf(args[10])
hist(log2(tpm_filter_mean$mean_tpm), breaks = seq(-25,20,0.5),main = tissue)
abline(v= 1,col = 2,lty=3,lwd =5)
text(x =-2.3, y = 800, labels = "x = 1" ,col = "red",cex = 1.5)
dev.off()

library(preprocessCore)
tpm_normquan <- normalize.quantiles(as.matrix(tpm_filter))
nqt<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
tpm_nqt<-t(apply(tpm_normquan,MARGIN=1,FUN=nqt))
colnames(tpm_nqt)<-colnames(tpm_filter)
rownames(tpm_nqt)<-rownames(tpm_filter)
save(tpm_nqt,file = args[5])

library(PCAForQTL)
expr<- t(tpm_nqt)
prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE)
PCs<-prcompResult$x 

importanceTable<-summary(prcompResult)$importance
PVEs<-importanceTable[2,]
sum(PVEs) #Theoretically, this should be 1.
pdf()
plot(PVEs,xlab="PC index",ylab="PVE")
dev.off()
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)

set.seed(2) #No need to change the RNG type since mc.cores will need to be 1.
resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05,
                              mc.cores=1)
print(resultRunBE$numOfPCsChosen)

K_elbow<-resultRunElbow 
K_BE<-resultRunBE$numOfPCsChosen 
K_GTEx<-60 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
#plot selected K
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(K_elbow,K_BE,K_GTEx),
                         titleText=tissue)
ggplot2::ggsave(args[6], width=16,height=11,unit="cm")
#know Covariate
eigenvec <- read.table(args[8])
eigenvec$V1 <- as.character(eigenvec$V1)
rownames(eigenvec) <- eigenvec$V1
cov1 <- eigenvec[as.character(group$WGS_ID),3:(n+2)]
colnames(cov1) <- c(paste0("PC",seq(1:n)))
unique(rownames(cov1)==group$WGS_ID)
rownames(cov1) <- rownames(expr)
gender=gsub("Female","1",group$gender)
gender=gsub("Male","2",gender)
knownCovariates <- cbind(cov1,gender)
knownCovariates$gender <- as.numeric(knownCovariates$gender)
identical(rownames(knownCovariates),rownames(expr)) #TRUE is good.

colnames(knownCovariates)[1:n]<-paste0("genotypePC",1:n)
PCsTop<-PCs[,1:K_BE] 
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,PCsTop,unadjustedR2_cutoff=0.9)
PCsTop<-scale(PCsTop) #Optional. Could be helpful for avoiding numerical inaccuracies.
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)
save(covariatesToUse,PCsTop,gender,K_BE,file=args[7])

gene.exp.nqt <- t(tpm_nqt)
gene.exp.adj <- gene.exp.nqt
for (i in 1:ncol(gene.exp.nqt)) {
  data = as.data.frame(cbind(gene.exp.nqt[, i], covariatesToUse))
  colnames(data) <- c("TPM", colnames(covariatesToUse))
  this.fit <- lm(TPM ~ ., data = data)
  gene.exp.adj[, i] <- coefficients(this.fit)[1] + residuals(this.fit)
}
save(gene.exp.adj, file = args[9])

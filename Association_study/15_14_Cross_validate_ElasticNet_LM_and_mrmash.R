
#==========================================================================================
#=====================Cross validate--Elastic Net Regression=================================================
#==========================================================================================
gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
Predict_output="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/Predixcan/Predict_output"
main_path="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate"

pdf(paste0(main_path,"/cor_predict_vs_true.pdf"))
for(tissue in c("Mu","Li","AF","BF")){
cor_result=NULL
for(times in 1:5){
test_expr <- read.table(paste0(gene_input,"/",tissue,"_tpm_nqt_test",times,".txt"),header=T,row.names=1,sep="\t")
Predict_expr <- read.table(paste0(Predict_output,"/",tissue,"_18295SNP_predict_cross",times,".txt"),header=T,sep="\t")
rownames(Predict_expr) <- Predict_expr[,1]
Predict_expr <- Predict_expr[,-c(1:2)]
all_gene=intersect(colnames(Predict_expr), colnames(test_expr))
for(gene in all_gene){
 test_expr$ID=rownames(test_expr)
 Predict_expr$ID=rownames(Predict_expr)
 test_expr_single <- test_expr[,c("ID",gene)]
 Predict_expr_single <- Predict_expr[,c("ID",gene)]
 merge_test_predict <- merge(test_expr_single,Predict_expr_single,by = "ID")
 cor_predict_vs_true=cor(merge_test_predict[,2],merge_test_predict[,3])
 cor_result_single <- data.frame(name=gene,predict_vs_true=cor_predict_vs_true)
 cor_result <- rbind(cor_result,cor_result_single)
cat(gene,times,"\n")
}}

#------------------------picture
aggr_mean <- aggregate(cor_result,by=list(cor_result$name),FUN=mean)
h <- hist(aggr_mean$predict_vs_true,plot = FALSE)
colors <- colorRampPalette(c("#98F5FF","#172869"))(length(h$counts))
mean=round(mean(na.omit(aggr_mean$predict_vs_true)),digits = 2)

par(mfrow = c(1,1),omi=c(0.5,0.5,0.5,0.5),mar=c(4,4,4,4),family = "sans")
hist(aggr_mean$predict_vs_true,labels = T,
     col=colors, family = "sans",
     cex.lab = 1.2, las = 1,cex.main=1,main=paste0(tissue,"--Elastic Net","\n","mean=",mean),
     ylim=c(0,max(h$counts)+20),
     xlab = "")
mtext("Pearson correlation coefficient",side = 1, line = 2, adj = 0.5,cex = 1.2)
box(lty = 1, lwd = 2,bty="l")
}
dev.off()

#-------------------示例
tissue="BF"
aggr_mean <- aggregate(cor_result,by=list(cor_result$name),FUN=mean)
h <- hist(aggr_mean$predict_vs_true,plot = FALSE)
colors <- colorRampPalette(c("#98F5FF","#172869"))(length(h$counts))
mean=round(mean(na.omit(aggr_mean$predict_vs_true)),digits = 2)

pdf("/mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Picture/Cross_validate_Elastic_Net_Example.pdf",width = 6,height = 6)
par(mfrow = c(1,1),omi=c(0.5,0.5,0.5,0.5),mar=c(4,4,4,4),family = "sans")
hist(aggr_mean$predict_vs_true,labels = T,
     col=colors, family = "sans",
     cex.lab = 1.2, las = 1,cex.main=1,main=paste0("Backfat adipose --Elastic Net","\n","mean=",mean),
     ylim=c(0,max(h$counts)+20),
     xlab = "")
mtext("Pearson correlation coefficient",side = 1, line = 2, adj = 0.5,cex = 1)
box(lty = 1, lwd = 2,bty="l")
dev.off()

#==========================================================================================
#=====================Cross validate--LM=================================================
#==========================================================================================
source("/mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh");

gene_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/gene_input"
snp_input="/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate/snp_input"
load(paste0(expression_group,"/Four_Tissues_group_rm_outliers.RData"))

sample_list <- read.table("/mnt/research/qgg/wu/YFJH/eQTL_mapping/expression_group/481_all_keep_indivi_WGSid.txt", header = F)

set.seed(3)
gene_sample <- sample(sample_list[,1])
k=floor(nrow(sample_list)/5)

for(times in 1:5){
if(times < 5){
gene_sample2 <- gene_sample[(k*(times-1)+1):(k*times)]
cat(k*(times-1)+1,k*times,"/n")
} else{
gene_sample2 <- gene_sample[(k*(times-1)+1):nrow(sample_list)]
cat(k*(times-1)+1,nrow(sample_list),"/n")
}

pdf(paste0("/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/PredictDB/cross_validate","/LM_cor_predict_vs_true.pdf"))

#expr of each tissue 
for(tissue in c("Mu","Li","AF","BF")){
load(paste0(PCAForQTL,"/",tissue,"_PCAForQTL_50PCs.RData"))
load(paste0(expression_group,"/",tissue,"_tpm.RData"))
#chr 1-18
gene_list <- read.table(paste0(eQTL_egene_summarize,"/",tissue,"_ms_uniq_cis_eGene.list")) # chr 1-18
#group
group=get(paste0(tissue,"_group"))

if (identical(as.character(group$Sample_id),colnames(tpm_nqt))){
  colnames(tpm_nqt) <- group$WGS_ID
} else {
  tpm_nqt <- tpm_nqt[,group$Sample_id]
  colnames(tpm_nqt) <- group$WGS_ID
}
tpm_nqt_filter <- tpm_nqt[gene_list$V1,] #cis-eGene(FRD=0.05(H2),Autosomes,cis-eQTL)
tpm_nqt_filter <- t(tpm_nqt_filter)

#get train set and text set 
tpm_nqt_test <- tpm_nqt_filter[rownames(tpm_nqt_filter) %in% gene_sample2, ]
tpm_nqt_filter2 <- tpm_nqt_filter[!(rownames(tpm_nqt_filter) %in% gene_sample2), ]

#-------------------------covariates
load(paste0(PCAForQTL,"/",tissue,"_PCAForQTL_50PCs.RData"))
if (identical(as.character(group$Sample_id),rownames(covariatesToUse))){
  rownames(covariatesToUse) <- group$WGS_ID
} else {
  covariatesToUse <- covariatesToUse[group$Sample_id,]
  rownames(covariatesToUse) <- group$WGS_ID
}

covariatesToUse <- t(covariatesToUse) 
covariatesToUse <- covariatesToUse[,!(colnames(covariatesToUse) %in% rownames(tpm_nqt_test))]
rownames(covariatesToUse) <- seq(1,nrow(covariatesToUse),1)

if(identical(colnames(covariatesToUse),rownames(tpm_nqt_filter2))){
write.table(tpm_nqt_filter2,file=paste0(gene_input,"/",tissue,"_transformed_expression",times,".txt"),col.names=T,row.names=T,sep="\t")
write.table(tpm_nqt_test,file=paste0(gene_input,"/",tissue,"_tpm_nqt_test",times,".txt"),col.names=T,row.names=T,sep="\t")
write.table(covariatesToUse,file=paste0(gene_input,"/",tissue,"_covariates",times,".txt"),col.names=T,row.names=T,sep="\t")
}
#-------------------------genotype
genotype <- read.table(paste0(snp_input,"/",tissue,"_genotype.txt"),header=T,check.names = F)
#genotype_train
genotype_train <- genotype[,c("varID",rownames(tpm_nqt_filter2))]
#genotype_text
genotype_text <- genotype[,c("varID",rownames(tpm_nqt_test))]

indivi_list <- data.frame(FID=rownames(tpm_nqt_test))
write.table(genotype_train,file=paste0(snp_input,"/",tissue,"_genotype_train",times,".txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(indivi_list,file=paste0(snp_input,"/",tissue,"_indivi_list_test",times,".txt"),col.names=F,row.names=F,sep=" ",quote=F)


#-------------------------lm predict 
#expr
adjusted.expr.file <- read.table(paste0(expression_group,"/",tissue,"_gene_exp_adj_PCAForQTL_FDR05.txt"), header = T, comment.char = "/", check.name = F) 
rownames(adjusted.expr.file) <- adjusted.expr.file$IID
adjusted.expr.filter <- adjusted.expr.file[,gene_list$V1] 
#get train set and text set 
adjusted.expr.train <- adjusted.expr.filter[!(rownames(adjusted.expr.filter) %in% gene_sample2), ]
adjusted.expr.test <- adjusted.expr.filter[rownames(adjusted.expr.filter) %in% gene_sample2, ]

#gentotype gene trio chr 1-18
gene_SNP_trio <- read.table(paste0(eQTL_egene_summarize,"/",tissue,"_eQTL_egene_summarize.txt"),header = T) # chr 1-18

#change SNP ID name  
rownames(genotype_train) <- substr(genotype_train$varID, 1, nchar(genotype_train$varID) - 8)
genotype_train2 <- genotype_train[,-1]
rownames(genotype_text) <- substr(genotype_text$varID, 1, nchar(genotype_text$varID) - 8)
genotype_text2 <- genotype_text[,-1]

cor_result=NULL
for (gene_line in 1:ncol(adjusted.expr.train)) {
gene_name=colnames(adjusted.expr.train)[gene_line]

#SNP
SNP_set <- gene_SNP_trio[gene_SNP_trio$gene %in% gene_name,]
#geno train
SNP_set_geno_train <- t(genotype_train2[rownames(genotype_train2) %in% SNP_set$snp,])
#geno test
SNP_set_geno_text <- as.data.frame(t(genotype_text2[rownames(genotype_text2) %in% SNP_set$snp,]))

#expr train
single.adjusted.expr.train <- as.data.frame(adjusted.expr.train[,gene_line])
rownames(single.adjusted.expr.train) <- rownames(adjusted.expr.train)
colnames(single.adjusted.expr.train) <- "expr"
#expr test
true_expr <- as.data.frame(adjusted.expr.test[,gene_name])
rownames(true_expr) <- rownames(adjusted.expr.test)
colnames(true_expr) <- "True_expr"
true_expr$ID <- rownames(true_expr)

if(identical(rownames(single.adjusted.expr.train), rownames(SNP_set_geno_train))){
data = as.data.frame(cbind(single.adjusted.expr.train, SNP_set_geno_train))
}else{
single.adjusted.expr.train <- single.adjusted.expr.train[rownames(SNP_set_geno_train),]
data = as.data.frame(cbind(single.adjusted.expr.train, SNP_set_geno_train))
}
#-----------fit
this.fit <- lm(expr ~ ., data = data)
predict_expr <- as.data.frame(predict(this.fit,newdata=SNP_set_geno_text))

colnames(predict_expr) <- "Predict_expr"
predict_expr$ID <- rownames(predict_expr)
expr_predict_vs_true <- merge(predict_expr, true_expr, by = "ID")
single_cor <- cor(expr_predict_vs_true$Predict_expr, expr_predict_vs_true$True_expr)
cor_result_single <- data.frame(name=gene_name,cis_SNP_num=length(SNP_set$snp),predict_vs_true=single_cor)

cor_result <- rbind(cor_result,cor_result_single)
}

h <- hist(cor_result$predict_vs_true,plot = FALSE)
colors <- colorRampPalette(c("#98F5FF","#172869"))(length(h$counts))
mean=round(mean(na.omit(cor_result$predict_vs_true)),digits = 2)

par(mfrow = c(1,1),omi=c(0.5,0.5,0.5,0.5),mar=c(4,4,4,4),family = "sans")
hist(cor_result$predict_vs_true,labels = T,
     col=colors, family = "sans",
     cex.lab = 1.2, las = 1,cex.main=1,main=paste0(tissue,"--LM","\n","mean=",mean),
     ylim=c(0,max(h$counts)+20),
     xlab = "")
mtext("Pearson correlation coefficient",side = 1, line = 2, adj = 0.5,cex = 1.2)
box(lty = 1, lwd = 2,bty="l")
dev.off()
}

#==========================================================================================
#=====================Cross validate--MR.mash=================================================
#==========================================================================================
pdf("/mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Picture/Cross_validate_mr.mash.pdf",width = 8,height = 6)
cor_raw <- read.table("/mnt/research/qgg/wu/YFJH/eQTL_mapping/TWAS/MR.mash/all_gene/output/mr.mash.alpha.predict_cor.txt",header = F)
colnames(cor_raw) <- c("cor_result_single","gene","fit_vs_true.cor","predict_vs_true.cor","cis_snp_count")
aggr_mean <- aggregate(cor_raw,by=list(cor_raw$gene),FUN=mean)
h <- hist(aggr_mean$predict_vs_true.cor,plot = FALSE)
colors <- colorRampPalette(c("#98F5FF","#172869"))(length(h$counts))
mean=round(mean(na.omit(aggr_mean$predict_vs_true.cor)),digits = 2)

par(mfrow = c(1,1),omi=c(0.5,0.5,0.5,0.5),mar=c(4,4,4,4),family = "sans")
hist(aggr_mean$predict_vs_true.cor,labels = T,
     col=colors, family = "sans",
     cex.lab = 1.2, las = 1,cex.main=1,main=paste0("mr.mash","\n","mean=",mean),
     ylim=c(0,max(h$counts)+20),
     xlab = "")
mtext("Pearson correlation coefficient",side = 1, line = 2, adj = 0.5,cex = 1.2)
box(lty = 1, lwd = 2,bty="l")

dev.off()

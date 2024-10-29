
args <- commandArgs(TRUE)
library(doMC)
start=args[1]
end=args[2]
registerDoMC(as.integer(args[3]))
source("/mnt/ufs18/rs-015/qgg/wu/YFJH/eQTL_mapping/code/Absolute_path.sh"); 
input=args[4]
cis_genotype=args[5]
cis_gene_list=args[6]
picture=args[7]
output=args[8]
#--------ciseQTL of common gene
common_gene_ciseQTL_list <- read.table(paste0(input,"/Four_tissues_common_gene_cis_SNP.ID"),header=F); 

#--------common gene list
common_gene_list <- read.table(paste0(expression_group,"/all_tissue_9614common_genename.txt"),header=F)
gene_anno <- read.table(paste0(expression_group,"/Sus_scrofa_11.1.104_gene_anno.txt"),header=T)
#gene_filter <- gene_anno[gene_anno$gene_id %in% common_gene_list$V1 & gene_anno$chr %in% c(1:18,"X","Y","MT"),] #9454 gene keeped;1-18:9104; X:316; Y:7; MT:27
gene_filter <- gene_anno[gene_anno$gene_id %in% common_gene_list$V1 & gene_anno$chr %in% c(1:18),] #9104 gene keeped;

#--------expr
load(paste0(expression_group,"/Four_Tissues_group_rm_outliers.RData"))
train_set_expr= read.table(paste0(expression_group,"/481_all_keep_indivi_WGSid.txt"), row.names=1)

#--------
mr.mash.alpha <- foreach(i = start:end, .combine = rbind) %dopar% {  

#for (gene_line in start:end){
gene_name=gene_filter[i,1]
#gene_name="ENSSSCG00000040140"
allsnp_genotype <- read.table(paste0(cis_genotype,"/",gene_name,"_train_dataset_cisSNP_for_mrmash.raw"),header=T,check.names=F); 
rownames(allsnp_genotype) <- allsnp_genotype$FID
common_gene_ciseQTL_genotype <- allsnp_genotype[,-c(1:6)]
geno.var <- apply(common_gene_ciseQTL_genotype, 2, var, na.rm = TRUE)
common_gene_ciseQTL_genotype <- common_gene_ciseQTL_genotype[, which(geno.var > 0.05)]

#--------------------------- train set expr
tissue=c("Mu","Li","AF","BF")
for(tissue_line in 1:4){
  load(paste0(expression_group,"/",tissue[tissue_line],"_gene_exp_adj_PCAForQTL_50PCs.RData"))
  eval(parse(text= paste0(tissue[tissue_line],"_gene_expr <- gene.exp.adj")))
  gene_expr <- get(paste0(tissue[tissue_line],"_gene_expr"))
  common_gene_expr <- as.data.frame(gene_expr[,colnames(gene_expr) %in% gene_name])
  colnames(common_gene_expr) <- tissue[tissue_line]
  common_gene_expr$ID <- rownames(gene_expr)
#group
  group <- get(paste0(tissue[tissue_line],"_group"))
  
  common_gene_expr_temp <- merge(group[,c("Sample_id","WGS_ID")],common_gene_expr,by.x="Sample_id",by.y="ID")
  rownames(common_gene_expr_temp) <- common_gene_expr_temp$WGS_ID
  common_gene_expr_single <- common_gene_expr_temp[,c(2,3)]
  train_set_expr <- merge(train_set_expr,common_gene_expr_single, by.x = "V2", by.y = "WGS_ID",all=T)
  #train_set_expr <- merge(train_set_expr,common_gene_expr_single, by.x = "V2", by.y = "WGS_ID")
 }
 rownames(train_set_expr)  <- train_set_expr[,1]
 train_set_expr <- train_set_expr[,-1]

 if(identical(rownames(train_set_expr), rownames(common_gene_ciseQTL_genotype))){
 dat_list <- list(X=common_gene_ciseQTL_genotype,Y=train_set_expr)
 }else{
common_gene_ciseQTL_genotype_sub <- common_gene_ciseQTL_genotype[rownames(train_set_expr),]
#test_set_common_gene_ciseQTL_genotype_sub <- test_set_common_gene_ciseQTL_genotype[,colnames(common_gene_ciseQTL_genotype)]
dat_list <- list(X=common_gene_ciseQTL_genotype_sub,Y=train_set_expr)
}
eval(parse(text= paste0(gene_name,"_gene_expr <- dat_list")))
save(dat_list,file=paste0(cis_gene_list,"/",gene_name,"_single_gene_list.RData"))

#-------------------------mr.mash.alpha
library(mr.mash.alpha)

#-------------------------cross validation
set.seed=1
gene_sample <- sample(1:nrow(dat_list$X))
k=ceiling(nrow(dat_list$X)/5)

cor_result=NULL
for(times in 1:5){
if(times < 5){
gene_sample2 <- gene_sample[(k*(times-1)+1):(k*times)]
cat(k*(times-1)+1,k*times,"/n")
} else{
gene_sample2 <- gene_sample[(k*(times-1)+1):nrow(dat_list$X)]
cat(k*(times-1)+1,nrow(dat_list$X),"/n")
}

#set.seed(1)

Ytrain <- as.matrix(dat_list$Y[-c(gene_sample2), ])
Xtrain <- as.matrix(dat_list$X[-c(gene_sample2), ])
Ytest <- as.matrix(dat_list$Y[c(gene_sample2), ])
Xtest <- as.matrix(dat_list$X[c(gene_sample2), ])

 
univ_sumstats <- compute_univariate_sumstats(Xtrain, Ytrain,
                   standardize=TRUE, standardize.response=FALSE)

grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2

#data-driven covariance matrices
S0_data_driven <- compute_data_driven_covs(sumstats=univ_sumstats,n_pcs=3)
S0 <- expand_covs(S0_data_driven, grid, zeromat=TRUE)
 
####Fit mr.mash
fit <- mr.mash(Xtrain, Ytrain, S0, update_V=TRUE, standardize = TRUE)

# Predict the multivariate outcomes in the test set using the fitted model.
Ytest_est <- predict(fit,Xtest)

#fit_list[]<-

cor_fit_vs_true <- cor.test(as.vector(fit$fitted),as.vector(Ytrain))$estimate
cor_predict_vs_true <- cor.test(as.vector(Ytest_est),as.vector(Ytest))$estimate
cis_snp_count <- ncol(Xtrain)

cor_result_single <- c(name=gene_name,fit_vs_true=cor_fit_vs_true,predict_vs_true=cor_predict_vs_true,cis_snp_count=cis_snp_count)
cor_result <- rbind(cor_result,cor_result_single)
cat(i,gene_name,times,"\n")
}
cor_result
}
write.table(mr.mash.alpha,file=paste0(output,"/",start,"_",end,"_mr.mash.alpha.predict_cor.txt"),quote=F)


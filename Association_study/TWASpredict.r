#==========================================================================================================================================
#===========================================================  Predict expression===========================================================
#==========================================================================================================================================

args=commandArgs(TRUE)
test_snp_input_file=args[1]
output_path=args[2]
suffix=args[3]

if(args[4]=="all"){
Tissue=c("Mu","Li","AF","BF")
}else{
Tissue=args[4]
}


input_path="/TWAS/LM/DLY_traits"

for(tissue in Tissue){
train_set_expr=paste0(input_path,"/",tissue,"_training_set_for_predict_expr.RData")
load(train_set_expr)

#test genotype
test_genotype <- read.table(test_snp_input_file,header=T,check.names = F)
test_genotype <- test_genotype[,-c(1,3:6)]
rownames(test_genotype) <- test_genotype[,1]
colnames(test_genotype) <- substr(colnames(test_genotype), 1, nchar(colnames(test_genotype)) - 4)
#colnames(test_genotype) <- substr(colnames(test_genotype), 1, nchar(colnames(test_genotype)) - 2)
#-------------------------start
all_predict_expr <- data.frame(IID=rownames(test_genotype))
for (gene_line in 1:length(train_set_list)) {
gene_name=names(train_set_list)[gene_line]
#train_set
train_set <- train_set_list[[gene_line]]
snp_set=colnames(train_set)[-1]
if(unique(is.na(match(snp_set,colnames(test_genotype))))){ 
predict_expr=NA
} else{
#extract genotype
column_indices2 <- match(snp_set,colnames(test_genotype)) 
column_indices2 <- column_indices2[!is.na(column_indices2)]

SNP_set_geno_text <- as.data.frame(test_genotype[,column_indices2])
colnames(SNP_set_geno_text) <- colnames(test_genotype)[column_indices2]
rownames(SNP_set_geno_text) <- rownames(test_genotype)
#
train_set <- as.data.frame(train_set[,c("expr",colnames(test_genotype)[column_indices2])])
#fit
gene_name_fit <- lm(expr ~ .,data=train_set)
#predict
predict_expr <- as.data.frame(predict(gene_name_fit,newdata=SNP_set_geno_text))
colnames(predict_expr) <- gene_name
predict_expr$IID <- rownames(predict_expr)

all_predict_expr <- merge(all_predict_expr, predict_expr, by = "IID")
}
}
write.table(all_predict_expr, file = paste0(output_path,"/",suffix,"_",tissue,"_predict_expr.txt"), row.names = F, col.names = T, quote = T)

}


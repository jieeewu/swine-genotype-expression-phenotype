===============================================================================================================
==================================== enrichment_of_eQTLs_in_QTL_summary =======================================
===============================================================================================================
args=commandArgs(TRUE)
pheno_file=args[1]
observe_file=args[2]
intera_path=args[3]
breed=args[4]
suffix=args[5]
output=args[6]

pheno_name=read.table(pheno_file)

observe <- read.table(observe_file)
colnames(observe) <- c("tissue","pheno","observe_count")

iter.mean.all=NULL
for (pheno in c(pheno_name$V1,"all_pheno")){
intera_single_pheno <- read.table(paste0(intera_path, "/", breed, "_", pheno,suffix))
colnames(intera_single_pheno) <- c("breed","pheno","tissue","times","intera_count")
iter.mean <- aggregate(intera_single_pheno[,5], by = list(intera_single_pheno$tissue), FUN = mean)
iter.mean$pheno <- pheno
iter.mean <- iter.mean[,c(3,1,2)]
colnames(iter.mean) <-  c("pheno","tissue","iter.mean.count")

max.count <- aggregate(intera_single_pheno, by = list(intera_single_pheno$tissue), FUN = max)[,c(3,4,6)]
colnames(max.count )[3]="max.count "
min.count <- aggregate(intera_single_pheno, by = list(intera_single_pheno$tissue), FUN = min)[,c(3,4,6)]
colnames(min.count)[3]="min.count"
max.min.count <- merge(max.count, min.count, by = c("pheno", "tissue"))
mean.max.min.count <- merge(iter.mean,max.min.count , by = c("pheno", "tissue"))
iter.mean.all <- rbind(iter.mean.all,mean.max.min.count)
}
reult_merge <- merge(observe,iter.mean.all,by = c("pheno", "tissue"))
reult_merge$mean_ratio <- reult_merge$observe_count / reult_merge$iter.mean.count
reult_merge$min.ratio <- reult_merge$observe_count / reult_merge$min.count
reult_merge$max.ratio <- reult_merge$observe_count / reult_merge$max.count
write.table(reult_merge,file=output,quote=F, row.names =F)




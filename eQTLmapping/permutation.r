# =================
# = Permutation =
# =================
args <- commandArgs(TRUE)
tissue <- args[1]
all_indivi <- read.table(args[3], header = F)
load(args[2])
set.seed(2)
all_indivi$V1 <- as.character(all_indivi$V1)
perm.seq <- matrix(ncol = 100, nrow = 481)
for (i in 1:100) { perm.seq[, i] <- sample(all_indivi$V1) }
load(args[5])
group <- get(paste0(tissue, "_group"))
unique(group$Sample_id==rownames(gene.exp.adj))
perm_gene <- matrix(ncol = 100, nrow = nrow(gene.exp.adj))
perm.seq2 <- apply(perm.seq, 2, as.character)
  
for (j in 1:100) {
  raw_ID <- as.character(group$WGS_ID)
  perm_gene[, j] <- raw_ID[na.omit(match(perm.seq2[, j], raw_ID))]
}
for (m in 1:ncol(gene.exp.adj)) {
  perm_pheno <- matrix(ncol = 100, nrow = nrow(gene.exp.adj))
  for (k in 1:100) { 
    perm_pheno[,k] <- gene.exp.adj[perm_gene[,k],m] 
  }
  perm_pheno2 <- cbind.data.frame(group$WGS_ID,group$WGS_ID, perm_pheno)
  colnames(perm_pheno2) <- c("#FID", "IID", paste0(colnames(gene.exp.adj)[m],"_",seq(1,100,1)))
  write.table(perm_pheno2, file = paste0(tissue, "/", colnames(gene.exp.adj)[m],".pheno"),quote = F, col.names = T,row.names = F,sep = " ")
  cat(m, "\n")
  }
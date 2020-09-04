## convert to normal format（for gcta or ldak GRM）
library(reshape2)
##for gcta
#tmp <- read.table(gzfile("snp.gcta.grm.gz"), header = F, stringsAsFactors = F)
#ids <- read.table("snp.gcta.grm.id", header = F, stringsAsFactors = F)
##for ldak
tmp <- read.table(gzfile("snp.ldak.weight.grm.gz"), header = F, stringsAsFactors = F)
ids <- read.table("snp.ldak.weight.grm.id", header = F, stringsAsFactors = F)
tmp <- tmp[,c(1,2,4)]
result_matrix <- acast(tmp, V1~V2, value.var="V4", drop = F)
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
result_full <- makeSymm(result_matrix)
#把对角线设为2（对角线值为1+近交系数的估计值）
#diag(result_full) <- 2
result_df <- as.data.frame(result_full)
row.names(result_df) <- ids$V2
colnames(result_df) <- ids$V2
#write.table(result_df, file = "gcta.kinship.txt", row.names = T, col.names = NA, sep = "\t", quote = F)
write.table(result_df, file = "ldak.weight.kinship.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

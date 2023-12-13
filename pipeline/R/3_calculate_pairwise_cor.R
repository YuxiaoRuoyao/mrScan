library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(purrr)
association_files <- unlist(snakemake@input[["association_files"]])
inst_files <- unlist(snakemake@input[["inst_files"]])
id_file <- read.csv(snakemake@input[["id_list"]])
out <- snakemake@output[["out"]]

id.list <- id_file$id
my_inst <- map(inst_files, function(f){
  readRDS(f)
})
beta_matrix <- map_dfc(association_files, function(f){
  readRDS(f)
})

names(my_inst) <- id.list
res <- data.frame(t(combn(id.list,2)))
res$cor <- map2(res$X1, res$X2, function(i, j) {
  sub_inst <- my_inst[c(i, j)]
  inst <- unique(unlist(sub_inst))
  sub_matrix <- beta_matrix %>% filter(SNP %in% inst)
  cor(sub_matrix[[i]],sub_matrix[[j]],use = "complete.obs")
}) %>% unlist()
res2 <- data.frame(X1 = res$X2, X2 = res$X1, cor = res$cor)
res_all <- rbind(res,res2)
df_matrix <- as.data.frame.matrix(xtabs(cor ~ ., res_all))
df_matrix <- abs(df_matrix)
saveRDS(list(df_pair = res, R_matrix = df_matrix),
        file = out)

library(dplyr)
library(purrr)
library(readr)
library(sumstatFactors)
library(stringr)

beta_files <- unlist(snakemake@input[["beta"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

ld <- purrr::map_dfr(1:22, function(c){
  read_table(ld_files[c])
})

M <- purrr:::map(1:22, function(c){
  read_lines(m_files[c])
}) %>% unlist() %>% as.numeric() %>% sum()

X <- map_dfr(beta_files, function(f){
  readRDS(f) %>%
    rename(SNP = snp) %>%
    inner_join(., ld)})
Z_hat <- X %>%
  select(ends_with(".z")) %>%
  as.matrix()
nms <- str_replace(colnames(Z_hat), ".z$", "")
SS <- X %>%
  select(ends_with(".ss")) %>%
  as.matrix()
N <- apply(SS, 2, median, na.rm=TRUE)
# Pairwise correlation
res <- data.frame(t(combn(nms,2)))
res$cor <- purrr::map2(res$X1, res$X2, function(i, j) {
  sub_Z <- Z_hat[,c(paste0(i,".z"), paste0(j,".z"))]
  sub_Z <- data.frame(sub_Z,check.names=F) %>% mutate(SNP = X$SNP)
  complete_sub_Z <- sub_Z[complete.cases(sub_Z[ ,-3]),]
  sub_SNP <- complete_sub_Z$SNP
  complete_sub_Z <- complete_sub_Z %>% select(-SNP) %>% as.matrix()
  sub_ld <- X %>% filter(SNP %in% sub_SNP) %>% pull(L2)
  sub_N <- N[c(paste0(i,".ss"),paste0(j,".ss"))]
  R <- R_ldsc(Z_hat = complete_sub_Z,
         ldscores = sub_ld,
         ld_size = M,
         N = sub_N,
         return_gencov = TRUE,
         return_cor = TRUE,
         make_well_conditioned = FALSE)
  R$Re[1,2]
}) %>% unlist()
res2 <- data.frame(X1 = res$X2, X2 = res$X1, cor = res$cor)
res_all <- rbind(res,res2)
df_matrix <- as.data.frame.matrix(xtabs(cor ~ ., res_all))
df_matrix <- abs(df_matrix)
saveRDS(list(df_pair = res, R_matrix = df_matrix),file=out)

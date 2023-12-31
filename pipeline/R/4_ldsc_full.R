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
R <- R_ldsc(Z_hat = Z_hat,
            ldscores = X$L2,
            ld_size = M,
            N = N,
            return_gencov = TRUE,
            return_cor = TRUE,
            make_well_conditioned = TRUE)
Re <- abs(R$Re)
Rg <- abs(R$Rg)
colnames(Re) <- rownames(Re) <- colnames(Rg) <- rownames(Rg) <- nms
saveRDS(list(Re = Re, Rg = Rg),file=out)

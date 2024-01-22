library(dplyr)
library(purrr)
library(ebnm)
library(esmr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
R <- readRDS(snakemake@input[["R"]])
R_type <- snakemake@params[["R_type"]]
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}

res <- MVMR_ESMR(beta_files = beta_files, R_matrix = R_matrix,
                 pval_threshold = pval_threshold)
saveRDS(res,file = out)
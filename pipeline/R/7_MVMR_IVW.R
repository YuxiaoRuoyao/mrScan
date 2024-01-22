library(dplyr)
library(TwoSampleMR)
library(purrr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
out <- snakemake@output[["out"]]

res <- MVMR_IVW(beta_files = beta_files,pval_threshold = pval_threshold)
saveRDS(res, file = out)

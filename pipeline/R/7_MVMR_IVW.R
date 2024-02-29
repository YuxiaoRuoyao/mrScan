library(dplyr)
library(TwoSampleMR)
library(purrr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
out <- snakemake@output[["out"]]

dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_IVW(dat = dat, pval_threshold = pval_threshold)
saveRDS(res, file = out)

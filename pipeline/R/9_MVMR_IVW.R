library(dplyr)
library(TwoSampleMR)
library(purrr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
out <- snakemake@output[["out"]]

dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_IVW(dat = dat, pval_threshold = pval_threshold, type = "local",
                effect_size_cutoff = effect_size_cutoff)
saveRDS(res, file = out)

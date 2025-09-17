library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(purrr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out <- snakemake@output[["out"]]

dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_IVW(dat = dat, pval_threshold = pval_threshold, type = "local",
                effect_size_cutoff = effect_size_cutoff,
                type_outcome = type_outcome, prevalence_outcome = prevalence_outcome)
saveRDS(res, file = out)


library(dplyr)
library(GRAPPLE)
library(purrr)
library(stringr)
suppressWarnings(library(mrScan))
library(ieugwasr)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
R <- readRDS(snakemake@input[["R"]])
R_type <- snakemake@params[["R_type"]]
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}
dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_GRAPPLE(dat = dat, R_matrix = R_matrix,
                    pval_threshold = pval_threshold, type = "local",
                    effect_size_cutoff = effect_size_cutoff,
                    type_outcome = type_outcome,prevalence_outcome = prevalence_outcome)
saveRDS(res,file = out)

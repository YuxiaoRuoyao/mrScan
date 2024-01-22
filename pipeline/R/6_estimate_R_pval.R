library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
p_thresh <- as.numeric(snakemake@params[["p_thresh"]])
out <- snakemake@output[["out"]]

Rpt <- estimate_R_pval(beta_files = beta_files,p_thresh=p_thresh)
saveRDS(Rpt, file=out)

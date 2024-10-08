library(dplyr)
library(purrr)
library(readr)
library(GFA)
library(stringr)
library(mrScan)
library(Matrix)

beta_files <- unlist(snakemake@input[["beta"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out <- snakemake@output[["out"]]

res <- ldsc_full(dat = beta_files, ld_files = ld_files, m_files = m_files)

saveRDS(res,file=out)

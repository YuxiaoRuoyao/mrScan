library(TwoSampleMR)
library(dplyr)
library(mrScan)
ex_dat1 <- readRDS(snakemake@input[["file1"]])
ex_dat2 <- readRDS(snakemake@input[["file2"]])
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
out <- snakemake@output[["out"]]

res <- bidirection_mr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, min_instruments = min_instruments,
                      effect_size_cutoff = effect_size_cutoff, R2_cutoff = R2_cutoff)
saveRDS(res,file = out)


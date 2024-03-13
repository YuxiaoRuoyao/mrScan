library(TwoSampleMR)
library(dplyr)
library(GRAPPLE)
library(mrScan)

ex_dat1 <- readRDS(snakemake@input[["file1"]]) # X/Y + M
ex_dat2 <- readRDS(snakemake@input[["file2"]]) # Z
ex_dat3 <- readRDS(snakemake@input[["file3"]]) # Z + M
ex_dat4 <- readRDS(snakemake@input[["file4"]]) # X/Y
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
out <- snakemake@output[["out"]]

res <- bidirection_mvmr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, ex_dat3 = ex_dat3,
                        ex_dat4 = ex_dat4, min_instruments = min_instruments,
                        effect_size_cutoff = effect_size_cutoff,
                        R2_cutoff = R2_cutoff,df_info = df_info)
saveRDS(res,file = out)

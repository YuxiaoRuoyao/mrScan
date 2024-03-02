library(TwoSampleMR)
library(dplyr)
library(GRAPPLE)
library(mrScan)

ex_dat1 <- readRDS(snakemake@input[["file1"]]) # X/Y + M
ex_dat2 <- readRDS(snakemake@input[["file2"]]) # Z
ex_dat3 <- readRDS(snakemake@input[["file3"]]) # Z + M
ex_dat4 <- readRDS(snakemake@input[["file4"]]) # X/Y
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
out <- snakemake@output[["out"]]

res <- bidirection_mvmr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, ex_dat3 = ex_dat3,
                        ex_dat4 = ex_dat4, min_instruments = min_instruments)
saveRDS(res,file = out)

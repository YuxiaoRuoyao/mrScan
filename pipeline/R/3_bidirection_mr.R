library(TwoSampleMR)
library(dplyr)
library(data.table)
library(mrScan)
ex_dat1 <- readRDS(snakemake@input[["file1"]])
ex_dat2 <- readRDS(snakemake@input[["file2"]])
method <- snakemake@params[["method"]]
over.dispersion <- as.logical(snakemake@params[["over_dispersion"]])
loss.function <- snakemake@params[["loss_function"]]
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
out <- snakemake@output[["out"]]

res <- bidirection_mr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, method = method,
                      over.dispersion = over.dispersion,loss.function = loss.function,
                      min_instruments = min_instruments)
saveRDS(res,file = out)


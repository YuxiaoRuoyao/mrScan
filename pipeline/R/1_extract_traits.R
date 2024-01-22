library(ieugwasr)
library(dplyr)
library(mrScan)

#source("R/helpers.R")
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
batch <- as.vector(snakemake@params[["batch"]])
pop <- snakemake@params[["population"]]
pval_x <- as.numeric(snakemake@params[["pval_instruments"]])
pval_z <- as.numeric(snakemake@params[["pval_traits"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
min_snps <- as.numeric(snakemake@params[["min_snps"]])
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
out <- snakemake@output[["out"]]

res <- extract_traits(id_exposure = id_exposure, id_outcome = id_outcome,
                      batch = batch, pop = pop, pval_x = pval_x,
                      pval_z = pval_z, r2 = r2, kb = kb, min_snps = min_snps,
                      min_instruments = min_instruments)
saveRDS(res,file = out)

library(ieugwasr)
library(dplyr)
library(reshape2)
library(mrScan)

#source("R/helpers.R")

res <- readRDS(snakemake@input[["file"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
method <- snakemake@params[["method"]]
res_cor <- readRDS(snakemake@input[["pairwise_cor"]])
extra_traits <- snakemake@params[["extra_traits"]]
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
Re <- res_cor$Re
Rg <- res_cor$Rg
df_pairs <- melt(Rg, value.name = "cor",varnames = c("id1","id2")) %>%
  filter(cor != 1)
df_matrix <- data.frame(Rg,check.names = FALSE)

res_unique <- unique_traits(id.list = id.list, df_info = df_info, R_matrix = df_matrix,
                            df_pairs = df_pairs, R2_cutoff = R2_cutoff, method = method,
                            extra_traits = extra_traits)
saveRDS(res_unique,file = out)

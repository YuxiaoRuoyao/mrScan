library(dplyr)
library(mrScan)

res <- readRDS(snakemake@input[["file"]])
res_bidirection <- readRDS(snakemake@input[["file_bidirection"]])
p_cutoff <- as.numeric(snakemake@params[["p_cutoff"]])
extra_traits <- snakemake@params[["extra_traits"]]
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
df_bidirection <- res_bidirection$df_bidirection

res_marginal <- marginal(id.list = id.list,df_info = df_info,
                         df_bidirection = df_bidirection,
                         p_cutoff = p_cutoff,extra_traits = extra_traits)
saveRDS(res_marginal,file = out)

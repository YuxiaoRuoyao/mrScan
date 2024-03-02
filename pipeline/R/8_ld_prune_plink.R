library(ieugwasr)
library(dplyr)
library(mrScan)

X <- readRDS(snakemake@input[["beta"]])
r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
clump_kb <- snakemake@params[["clump_kb"]]
ref_path  <- snakemake@params[["ref_path"]]
type <- snakemake@params[["ld_prioritization"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
out <- snakemake@output[["out"]]


res_X <- ld_prune_plink(X = X,r2_thresh=r2_thresh,clump_kb=clump_kb,
                        ref_path = ref_path,type=type,pthresh=pthresh)
saveRDS(res_X, file=out)

library(dplyr)
library(MRBEE)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
pleio_threshold <- as.numeric(snakemake@params[["pleio_p_thresh"]])
R <- readRDS(snakemake@input[["R"]])
R_type <- snakemake@params[["R_type"]]
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}
dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_MRBEE(dat = dat, R_matrix =  R_matrix,
                  pval_threshold = pval_threshold,
                  pleio_threshold = pleio_threshold,type = "local",
                  effect_size_cutoff = effect_size_cutoff)
saveRDS(res,file = out)

library(dplyr)
library(mrScan)
library(MVMR)

beta_files <- unlist(snakemake@input[["beta"]])
R <- readRDS(snakemake@input[["R"]])
df_info <- readRDS(snakemake@input[["file"]])$trait.info
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
R_type <- snakemake@params[["R_type"]]
extra_traits <- snakemake@params[["extra_traits"]]
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}

res <- strength_filter(beta_files = beta_files,R_matrix = R_matrix,df_info = df_info,
                       pval_threshold = pval_threshold, F_threshold = F_threshold,
                       extra_traits = extra_traits)
saveRDS(res,file = out)

library(dplyr)
library(mrScan)
library(MVMR)

beta_files <- unlist(snakemake@input[["beta"]])
R <- readRDS(snakemake@input[["R"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
R_type <- snakemake@params[["R_type"]]
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}
dat <- purrr::map_dfr(beta_files, readRDS)
res <- strength_filter(dat = dat,R_matrix = R_matrix,
                       df_info = df_info,
                       pval_threshold = pval_threshold, Filter = FALSE)
write.csv(res,file = out,row.names = F)

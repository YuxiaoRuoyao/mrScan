library(hdme)
library(dplyr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
id_exposure <- snakemake@params[["id_exposure"]]
df_info <- read.csv(snakemake@input[["trait_info"]])
radius_type <- snakemake@params[["radius_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
maxits <- as.numeric(snakemake@params[["maxits"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

dat <- purrr::map_dfr(beta_files, readRDS)
res_Lasso <- corrected_Lasso(dat = dat, type = "local",
                             id_exposure = id_exposure,df_info = df_info,
                             radius_type = radius_type,
                             seed = seed, maxits = maxits,
                             type_outcome = type_outcome,
                             prevalence_outcome = prevalence_outcome)
write.csv(data.frame(id = res_Lasso$id.list),file = out_id_list,row.names = F)
write.csv(res_Lasso$trait.info,file = out_trait_info,row.names = F)

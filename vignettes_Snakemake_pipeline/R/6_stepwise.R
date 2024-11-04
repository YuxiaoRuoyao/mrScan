library(stringr)
library(mrScan)

beta_files <- unlist(snakemake@input[["beta"]])
id_exposure <- snakemake@params[["id_exposure"]]
df_info <- read.csv(snakemake@input[["trait_info"]])
method <- snakemake@params[["method"]]
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

dat <- purrr::map_dfr(beta_files, readRDS)
res_stepwise <- stepwise(dat = dat,type = "local",
                         id_exposure = id_exposure, df_info = df_info,
                         type_outcome = type_outcome,
                         prevalence_outcome = prevalence_outcome,
                         method = method)
write.csv(data.frame(id = res_stepwise$id.list),file = out_id_list,row.names = F)
write.csv(res_stepwise$trait.info,file = out_trait_info,row.names = F)

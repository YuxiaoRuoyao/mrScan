library(dplyr)
library(mrScan)
library(MVMR)

beta_files <- unlist(snakemake@input[["beta"]])
R <- readRDS(snakemake@input[["R"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
R_type <- snakemake@params[["R_type"]]
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
extra_traits <- snakemake@params[["extra_traits"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]
out_df_strength <- snakemake@output[["out_df_strength"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}
dat <- purrr::map_dfr(beta_files, readRDS)
res <- strength_filter(dat = dat,dat_type = "local",R_matrix = R_matrix,
                       df_info = df_info,
                       pval_threshold = pval_threshold,
                       F_threshold = F_threshold,
                       effect_size_cutoff = effect_size_cutoff,
                       extra_traits = extra_traits, Filter = TRUE,
                       type_outcome = type_outcome, prevalence_outcome = prevalence_outcome)

write.csv(data.frame(id = res$id.list),file = out_id_list,row.names = F)
write.csv(res$trait.info,file = out_trait_info,row.names = F)
write.csv(res$df_strength,file = out_df_strength,row.names = F)

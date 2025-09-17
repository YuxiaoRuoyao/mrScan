library(dplyr)
library(MRBEE)
library(mrScan)
library(ieugwasr)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
pleio_threshold <- as.numeric(snakemake@params[["pleio_p_thresh"]])
R <- readRDS(snakemake@input[["R"]])
R_type <- snakemake@params[["R_type"]]
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}
df_info_exposure_outcome <- read.csv("/nfs/turbo/sph-jvmorr/CRP_project/pipeline/ebi-a-GCST90475667/results/df_info_exposure_outcome.csv")
ncase_outcome <- df_info_exposure_outcome %>% filter(id == "ebi-a-GCST90475667") %>% pull(ncase)
ncontrol_outcome <- df_info_exposure_outcome %>% filter(id == "ebi-a-GCST90475667") %>% pull(ncontrol)
dat <- purrr::map_dfr(beta_files, readRDS)
res <- MVMR_MRBEE(dat = dat, R_matrix =  R_matrix,
                  pval_threshold = pval_threshold,
                  pleio_threshold = pleio_threshold,type = "local",
                  effect_size_cutoff = effect_size_cutoff,
                  type_outcome = type_outcome, prevalence_outcome = prevalence_outcome)
saveRDS(res,file = out)

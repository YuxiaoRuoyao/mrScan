library(ieugwasr)
library(dplyr)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
batch <- as.vector(snakemake@params[["batch"]])
pop <- snakemake@params[["population"]]
pval_x <- as.numeric(snakemake@params[["pval_instruments"]])
pval_z <- as.numeric(snakemake@params[["pval_traits"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
min_snps <- as.numeric(snakemake@params[["min_snps"]])
type_exposure <- snakemake@params[["type_exposure"]]
type_candidate_traits <- snakemake@params[["type_candidate_traits"]]
file_path <- snakemake@params[["file_path"]]
ref_path <- snakemake@params[["ref_path"]]
df_candidate_traits <- snakemake@params[["df_candidate_traits"]]
out <- snakemake@output[["out"]]

if(type_exposure == "local" & is.na(file_path)){
  stop("Please provide local path of the exposure data in config file!")
}else if(type_exposure == "local" & is.na(ref_path)){
  stop("Please provide LD reference path in config file!")
}else if(type_candidate_traits == "local" & is.na(df_candidate_traits)){
  stop("Please provide info file of candidate traits in config file!")
}

if(type_candidate_traits == "local" & (!is.na(df_candidate_traits))){
  df_info <- read.csv(df_candidate_traits,header = TRUE)
  file_list <- df_info$path
  trait_list <- df_info$trait_ID
  snp_name_list <- df_info$snp
  beta_hat_name_list <- df_info$beta_hat
  se_name_list <- df_info$se
  p_value_name_list <- df_info$p_value
}
res <- extract_traits(id_exposure = id_exposure, pval_x = pval_x, pval_z = pval_z,
                      pop = pop, batch = batch,
                      r2 = r2, kb = kb,
                      min_snps = min_snps,
                      type_exposure = type_exposure,
                      type_candidate_traits = type_candidate_traits,
                      file_path = file_path, ref_path = ref_path,
                      file_list = file_list,trait_list = trait_list,
                      snp_name_list = snp_name_list,
                      beta_hat_name_list = beta_hat_name_list,
                      se_name_list = se_name_list,
                      p_value_name_list = p_value_name_list)
saveRDS(res,file = out)

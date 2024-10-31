library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
library(mrScan)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ5dXhpYW93QHVtaWNoLmVkdSIsImlhdCI6MTcyOTY0MDQ5NSwiZXhwIjoxNzMwODUwMDk1fQ.YZA2h58w7WgJjxWzNnKKnYRz84c7Gf40FVTP4aKMdYpa84DkcPk-IoHFeJgjMEeIhwGEC7gXnrtb6MTW0pp24jEESZFVX7_hIg5k0UnX3q6Y4o9GhjppNo5dZFof9ZRZmMzu43WCrUj_nYNcWL-k7vGaXhpGaS2AC1sAk24RLDHY23HYyPEQ76O0Cs4hBDjF2veNj8XAbV_6ThSLX8NA2g4ygd3FiZjlXk8cEPeEF49thNfNQWOJD3rcfI52V5AIn1VJk5gIN4PD-fm1GLWHou02EDdlVt_-sewZlPjmkpu7jVMjP9OZnppwIYj4w6Bi6uf4uYeFzMyiIrP6bKYPng")

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
print(batch)
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

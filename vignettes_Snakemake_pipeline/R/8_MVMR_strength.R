library(dplyr)
library(mrScan)
library(MVMR)
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ5dXhpYW93QHVtaWNoLmVkdSIsImlhdCI6MTcyOTY0MDQ5NSwiZXhwIjoxNzMwODUwMDk1fQ.YZA2h58w7WgJjxWzNnKKnYRz84c7Gf40FVTP4aKMdYpa84DkcPk-IoHFeJgjMEeIhwGEC7gXnrtb6MTW0pp24jEESZFVX7_hIg5k0UnX3q6Y4o9GhjppNo5dZFof9ZRZmMzu43WCrUj_nYNcWL-k7vGaXhpGaS2AC1sAk24RLDHY23HYyPEQ76O0Cs4hBDjF2veNj8XAbV_6ThSLX8NA2g4ygd3FiZjlXk8cEPeEF49thNfNQWOJD3rcfI52V5AIn1VJk5gIN4PD-fm1GLWHou02EDdlVt_-sewZlPjmkpu7jVMjP9OZnppwIYj4w6Bi6uf4uYeFzMyiIrP6bKYPng")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
beta_files <- unlist(snakemake@input[["beta"]])
R <- readRDS(snakemake@input[["R"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
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
dat <- purrr::map_dfr(beta_files, readRDS)
res <- strength_filter(dat = dat,dat_type = "local",R_matrix = R_matrix,
                       df_info = df_info,
                       pval_threshold = pval_threshold, Filter = FALSE,
                       effect_size_cutoff = effect_size_cutoff,
                       type_outcome = type_outcome, prevalence_outcome = prevalence_outcome)
write.csv(res,file = out,row.names = F)

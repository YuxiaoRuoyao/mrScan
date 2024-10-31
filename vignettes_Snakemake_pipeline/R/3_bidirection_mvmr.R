library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(GRAPPLE)
library(mrScan)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJ5dXhpYW93QHVtaWNoLmVkdSIsImlhdCI6MTcyOTY0MDQ5NSwiZXhwIjoxNzMwODUwMDk1fQ.YZA2h58w7WgJjxWzNnKKnYRz84c7Gf40FVTP4aKMdYpa84DkcPk-IoHFeJgjMEeIhwGEC7gXnrtb6MTW0pp24jEESZFVX7_hIg5k0UnX3q6Y4o9GhjppNo5dZFof9ZRZmMzu43WCrUj_nYNcWL-k7vGaXhpGaS2AC1sAk24RLDHY23HYyPEQ76O0Cs4hBDjF2veNj8XAbV_6ThSLX8NA2g4ygd3FiZjlXk8cEPeEF49thNfNQWOJD3rcfI52V5AIn1VJk5gIN4PD-fm1GLWHou02EDdlVt_-sewZlPjmkpu7jVMjP9OZnppwIYj4w6Bi6uf4uYeFzMyiIrP6bKYPng")

ex_dat1 <- readRDS(snakemake@input[["file1"]]) # X/Y + M
ex_dat2 <- readRDS(snakemake@input[["file2"]]) # Z
ex_dat3 <- readRDS(snakemake@input[["file3"]]) # Z + M
ex_dat4 <- readRDS(snakemake@input[["file4"]]) # X/Y
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
prevalence <- as.numeric(snakemake@params[["prevalence"]])
id_outcome <- snakemake@params[["id_outcome"]]
type_outcome <- snakemake@params[["type_outcome"]]
out <- snakemake@output[["out"]]

if(unique(ex_dat4$id.exposure) == id_outcome & type_outcome == "binary"){
  type_list <- c("binary","continuous")
  prevalence_list <- c(prevalence, NA)
}else{
  type_list <- c("continuous","continuous")
  prevalence_list <- c(NA, NA)
}

res <- bidirection_mvmr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, ex_dat3 = ex_dat3,
                        ex_dat4 = ex_dat4, min_instruments = min_instruments,
                        effect_size_cutoff = effect_size_cutoff,
                        R2_cutoff = R2_cutoff,df_info = df_info,
                        type_list = type_list, prevalence_list = prevalence_list)
saveRDS(res,file = out)

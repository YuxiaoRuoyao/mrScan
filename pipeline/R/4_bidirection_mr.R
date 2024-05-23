library(TwoSampleMR)
library(dplyr)
library(GRAPPLE)
library(mrScan)

ex_dat1 <- readRDS(snakemake@input[["file1"]])
ex_dat2 <- readRDS(snakemake@input[["file2"]])
min_instruments <- as.numeric(snakemake@params[["min_instruments"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
prevalence <- as.numeric(snakemake@params[["prevalence"]])
id_outcome <- snakemake@params[["id_outcome"]]
type_outcome <- snakemake@params[["type_outcome"]]
out <- snakemake@output[["out"]]

if(unique(ex_dat1$id.exposure) == id_outcome & type_outcome == "binary"){
  type_list <- c("binary","continuous")
  prevalence_list <- c(prevalence, NA)
}else if(unique(ex_dat2$id.exposure) == id_outcome & type_outcome == "binary"){
  type_list <- c("continuous","binary")
  prevalence_list <- c(NA,prevalence)
}else{
  type_list <- c("continuous","continuous")
  prevalence_list <- c(NA, NA)
}
res <- bidirection_mr(ex_dat1 = ex_dat1, ex_dat2 = ex_dat2, min_instruments = min_instruments,
                      effect_size_cutoff = effect_size_cutoff, R2_cutoff = R2_cutoff,
                      type_list = type_list,prevalence_list = prevalence_list)
saveRDS(res,file = out)


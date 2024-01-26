library(ieugwasr)
library(dplyr)
library(mrScan)

#source("R/helpers.R")
id_exposure <- snakemake@params[["id_exposure"]]
batch <- as.vector(snakemake@params[["batch"]])
pop <- snakemake@params[["population"]]
pval_x <- as.numeric(snakemake@params[["pval_instruments"]])
pval_z <- as.numeric(snakemake@params[["pval_traits"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
min_snps <- as.numeric(snakemake@params[["min_snps"]])
out <- snakemake@output[["out"]]

res <- extract_traits(id_exposure = id_exposure, pval_x = pval_x, pval_z = pval_z,
                      pop = pop, batch = batch,
                      r2 = r2, kb = kb,
                      min_snps = min_snps,
                      type_exposure = "IEU",
                      type_candidate_traits = "IEU",
                      file_path = NULL, ref_path = NULL,file_list = NULL,
                      trait_list = NULL, snp_name_list = NULL,
                      beta_hat_name_list = NULL, se_name_list = NULL,
                      p_value_name_list = NULL)
saveRDS(res,file = out)

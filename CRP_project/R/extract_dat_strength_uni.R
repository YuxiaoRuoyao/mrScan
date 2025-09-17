library(TwoSampleMR)
library(ieugwasr)

ex_dat <- readRDS(snakemake@input[["ex_dat"]])
x1 <- snakemake@params[["trait1"]]
x2 <- snakemake@params[["trait2"]]
id_outcome <- snakemake@params[["id_outcome"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out <- snakemake@output[["out"]]

#ex_dat <- mv_extract_exposures(id_exposure = c(x1,x2), pval_threshold = inst_pval)
out_dat <- extract_outcome_data(ex_dat$SNP,id_outcome)
dat <- mv_harmonise_data(ex_dat,out_dat)
saveRDS(dat,file = out)

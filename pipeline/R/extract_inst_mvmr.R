library(TwoSampleMR)
x1 <- snakemake@params[["trait1"]]
x2 <- snakemake@params[["trait2"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out <- snakemake@output[["out"]]

df_inst <- mv_extract_exposures(id_exposure = c(x1,x2), pval_threshold = inst_pval)
saveRDS(df_inst,file = out)

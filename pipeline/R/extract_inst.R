library(TwoSampleMR)
x <- snakemake@params[["trait"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out <- snakemake@output[["out"]]

df_inst <- extract_instruments(outcomes = x, p1 = inst_pval)
saveRDS(df_inst,file = out)

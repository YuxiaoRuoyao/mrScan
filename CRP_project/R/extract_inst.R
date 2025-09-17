library(TwoSampleMR)
library(ieugwasr)

x <- snakemake@params[["trait"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out <- snakemake@output[["out"]]

df_inst <- tryCatch({
  extract_instruments(outcomes = x, p1 = inst_pval)
}, error = function(e) {
  NULL
})
saveRDS(df_inst,file = out)

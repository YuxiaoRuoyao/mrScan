library(ieugwasr)
library(dplyr)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

df <- gwasinfo(c(id_exposure,id_outcome))
write.csv(df,file = out, row.names=F)

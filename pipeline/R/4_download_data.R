library(stringr)
library(dplyr)
library(mrScan)
res <- readRDS(snakemake@input[["file"]])
df_harmonise <- read.csv(snakemake@input[["df_harmonise"]],header = F)
data_path <- snakemake@params[["path"]]
path_checkpoint <- snakemake@params[["checkpoint"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

id_list <- unique(c(res$id.list,id_exposure,id_outcome))
f <- download_gwas(id_list = id_list, df_harmonise = df_harmonise)
write.table(f,file = out,row.names = FALSE,
            col.names = FALSE, quote = FALSE)


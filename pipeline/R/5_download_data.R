library(stringr)
library(dplyr)
library(mrScan)
res <- read.csv(snakemake@input[["file"]])
df_harmonise <- read.csv(snakemake@input[["df_harmonise"]],header = F)
data_path <- snakemake@params[["path"]]
path_checkpoint <- snakemake@params[["checkpoint"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

id_list <- unique(c(res$id,id_exposure,id_outcome))
f <- download_gwas(id_list = id_list, df_harmonise = df_harmonise,
                   data_path = data_path, path_checkpoint = path_checkpoint)
write.table(f,file = out,row.names = FALSE,
            col.names = FALSE, quote = FALSE)


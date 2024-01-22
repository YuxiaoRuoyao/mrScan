library(stringr)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
res <- readRDS(snakemake@input[["file"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
method <- snakemake@params[["method"]]
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
mvdat_y <- mvdat$mvdat_y

res_stepwise <- stepwise(id_exposure = id_exposure,id.list = id.list,
                         df_info = df_info, mvdat_y = mvdat_y,
                         method = method)
saveRDS(res_stepwise,file = out)

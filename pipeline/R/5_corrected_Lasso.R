library(hdme)
library(dplyr)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
res <- readRDS(snakemake@input[["file"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
radius_type <- snakemake@params[["radius_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
maxits <- as.numeric(snakemake@params[["maxits"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
mvdat_y <- mvdat$mvdat_y

res_Lasso <- corrected_Lasso(id_exposure = id_exposure,id.list = id.list,df_info = df_info,
                             mvdat_y = mvdat_y,radius_type=radius_type,
                             seed = seed,maxits = maxits)
saveRDS(res_Lasso,file = out)

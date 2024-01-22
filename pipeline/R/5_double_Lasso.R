library(glmnet)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
res <- readRDS(snakemake@input[["file"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
lambda_type <- snakemake@params[["lambda_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
mvdat_x <- mvdat$mvdat_x
mvdat_y <- mvdat$mvdat_y

res_Lasso <- double_Lasso(id_exposure = id_exposure,id.list = id.list, df_info = df_info,
                          mvdat_x = mvdat_x, mvdat_y = mvdat_y,
                          lambda_type = "min",seed = 1)
saveRDS(res_Lasso,file = out)

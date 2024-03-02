library(hdme)
library(dplyr)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
id.list <- read.csv(snakemake@input[["id_list"]])$id
df_info <- read.csv(snakemake@input[["trait_info"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
radius_type <- snakemake@params[["radius_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
maxits <- as.numeric(snakemake@params[["maxits"]])
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

mvdat_y <- mvdat$mvdat_y

res_Lasso <- corrected_Lasso(id_exposure = id_exposure,id.list = id.list,df_info = df_info,
                             mvdat_y = mvdat_y,radius_type=radius_type,
                             seed = seed,maxits = maxits)

write.csv(data.frame(id = res_Lasso$id.list),file = out_id_list,row.names = F)
write.csv(res_Lasso$trait.info,file = out_trait_info,row.names = F)

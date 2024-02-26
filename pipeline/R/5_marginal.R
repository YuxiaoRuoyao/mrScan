library(dplyr)
library(mrScan)

id.list <- read.csv(snakemake@input[["id_list"]])$id
df_info <- read.csv(snakemake@input[["trait_info"]])
df_bidirection <- read.csv(snakemake@input[["file_bidirection"]])
p_cutoff <- as.numeric(snakemake@params[["p_cutoff"]])
extra_traits <- snakemake@params[["extra_traits"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

res_marginal <- marginal(id.list = id.list,df_info = df_info,
                         df_bidirection = df_bidirection,
                         p_cutoff = p_cutoff,extra_traits = extra_traits)

write.csv(data.frame(id = res_marginal$id.list),file = out_id_list,row.names = F)
write.csv(res_marginal$trait.info,file = out_trait_info,row.names = F)

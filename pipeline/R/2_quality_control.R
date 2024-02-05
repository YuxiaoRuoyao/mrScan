library(ieugwasr)
library(dplyr)
library(mrScan)
res_initial <- snakemake@input[["file"]]
nsnp_cutoff <- as.numeric(snakemake@params[["nsnp_cutoff"]])
pop <- snakemake@params[["population"]]
sex <- snakemake@params[["sex"]]
out_id_list <- snakemake@output[["id_list"]]
out_trait_info <- snakemake@output[["trait_info"]]

dat <- readRDS(res_initial)$trait.info

res <- quality_control(dat = dat,nsnp_cutoff = nsnp_cutoff,pop = pop,sex = sex)

write.csv(data.frame(id = res$id.list),file = out_id_list,row.names = F)
write.csv(res$trait.info,file = out_trait_info,row.names = F)

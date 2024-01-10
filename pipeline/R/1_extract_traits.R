library(ieugwasr)
library(dplyr)
#library(mrScan)

source("R/helpers.R")
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
batch <- as.vector(snakemake@params[["batch"]])
pop <- snakemake@params[["population"]]
pval_x <- as.numeric(snakemake@params[["pval_instruments"]])
pval_z <- as.numeric(snakemake@params[["pval_traits"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
min_snps <- as.numeric(snakemake@params[["min_snps"]])
out <- snakemake@output[["out"]]

batch1 <- c("ieu-a","ieu-b","ukb-b")
batch2 <- batch[!batch %in% batch1]
phe1 <- retrieve_traits(id_exposure, pval_x,pval_z,pop = pop, batch=batch1,
                       r2 = r2, kb = kb,
                       access_token = ieugwasr::check_access_token(),
                       min_snps =min_snps)
if(length(batch2) != 0){
    phe2 <- retrieve_traits(id_exposure, pval_x,pval_z,pop = pop, batch=batch2,
                       r2 = r2, kb = kb,
                       access_token = ieugwasr::check_access_token(),
                       min_snps =min_snps)
    id.list <- unique(c(phe1$phe$id,phe2$phe$id)) 
}else{
    id.list <- unique(phe1$phe$id)
}
# Delete X
id.list.initial <- id.list[!id.list %in% id_exposure]
df_trait <- gwasinfo(id.list.initial)
df_trait <- df_trait[,c("id","trait","sex","consortium","nsnp","note","sample_size",
                        "population","year")]
df_trait['status'] <- 'Initial List'
saveRDS(list(id.list=id.list.initial,trait.info=df_trait),
        file = out)


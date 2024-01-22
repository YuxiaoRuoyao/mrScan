library(ieugwasr)
library(dplyr)
library(data.table)
library(purrr)
library(mrScan)
mr_files <- unlist(snakemake@input[["mr_files"]])
id_file <- read.csv(snakemake@input[["id_list"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
id_exposure <- snakemake@params[["id_exposure"]]
sig_level <- as.numeric(snakemake@params[["sig_level"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
extra_trait <- snakemake@params[["extra_trait"]]
out <- snakemake@output[["out"]]

id.list <- id_file$id
res_mr <- map(mr_files, function(f){
  readRDS(f)
})
res <- do.call(Map, c(f = rbind, res_mr))
res_cor <- res$cor

res_downstream <- downstream_filter(id_exposure = id_exposure,id.list = id.list,
                                    df_info = df_info,res = res, sig_level = sig_level)
select_trait <- res_downstream$id.list
df_info <- res_downstream$trait.info
# delete high correlation traits with either X and Y
res_high_cor_XY <- filter_high_cor_XY(id_list = select_trait, df_info = df_info,
                   res_cor = res_cor,id_exposure = id_exposure,R2_cutoff = R2_cutoff)
select_trait <- res_high_cor_XY$id.list
df_info <- res_high_cor_XY$trait.info
if(extra_trait != "None"){
  select_trait <- c(select_trait,extra_trait)
  df_info[df_info$id %in% extra_trait,"status"] <- "select after downstream filtering"
}
saveRDS(list(id.list=select_trait,trait.info=df_info,df_bidirection = res_downstream$df_bidirection),
        file = out)

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
method <- snakemake@params[["method"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]
out_df_bidirection <- snakemake@output[["out_df_bidirection"]]

id.list <- id_file$id
empty_files <- c()
high_corr_files <- c()
for (f in mr_files) {
  rds_content <- readRDS(f)
  if (is_empty(rds_content)) {
    empty_files <- c(empty_files, f)
  } else if (is.null(rds_content$mr12) && is.null(rds_content$mr21)) {
    high_corr_files <- c(high_corr_files, f)
  }
}

empty_id <- empty_files %>% strsplit("_") %>% sapply(tail, 1) %>%
  strsplit(".RDS") %>% sapply(tail, 1) %>% unique()
high_corr_id <- high_corr_files %>% strsplit("_") %>% sapply(tail, 1) %>%
  strsplit(".RDS") %>% sapply(tail, 1) %>% unique()
df_info[df_info$id %in% empty_id,"status"] <- "delete due to not enough instruments"
df_info[df_info$id %in% high_corr_id,"status"] <- "delete since high cor with X or Y"
id.list <- id.list[!id.list %in% c(empty_id,high_corr_id)]


mr_files <- mr_files[!mr_files %in% c(empty_files,high_corr_files)]
res_mr <- map(mr_files, function(f){
  readRDS(f)
})
res <- do.call(Map, c(f = rbind, res_mr))

res_downstream <- downstream_filter(id_exposure = id_exposure,id.list = id.list,
                                    df_info = df_info,res = res, sig_level = sig_level,
                                    method = method)
select_trait <- res_downstream$id.list
df_info <- res_downstream$trait.info
if(extra_trait != "None"){
  select_trait <- c(select_trait,extra_trait)
  df_info[df_info$id %in% extra_trait,"status"] <- "select after downstream filtering"
}
write.csv(data.frame(id = select_trait),file = out_id_list,row.names = F)
write.csv(df_info,file = out_trait_info,row.names = F)
write.csv(res_downstream$df_bidirection, file = out_df_bidirection,row.names = F)

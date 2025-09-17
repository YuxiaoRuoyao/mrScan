library(dplyr)
library(stringr)
library(mrScan)
library(MVMR)

beta_files <- unlist(snakemake@input[["beta"]])
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
R <- readRDS(snakemake@input[["R"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
extra_traits <- unlist(snakemake@params[["extra_traits"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
out_flag <- snakemake@output[["out_flag"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_df_strength <- snakemake@output[["out_df_strength"]]

dat <- purrr::map_dfr(beta_files, readRDS)
beta_hat <- dat %>% select(ends_with(".beta"))
nms <- str_replace(names(beta_hat), ".beta", "")
id_list <- unique(c(id_outcome,id_exposure,nms))
R_matrix <- as.matrix(R$Re)
R_matrix <- R_matrix[id_list,id_list]

prevalence_exposure = rep(NA,length(id_list)-1)
if("ieu-a-31" %in% id_list){
  prevalence_exposure[which(id_list %in% "ieu-a-31")-1] = 0.000843
}
if("ieu-a-30" %in% id_list){
  prevalence_exposure[which(id_list %in% "ieu-a-30")-1] = 0.0003372 # CD 40%
}
if("ieu-a-32" %in% id_list){
  prevalence_exposure[which(id_list %in% "ieu-a-32")-1] = 0.0005058 # UC 60%
}
res <- strength_filter(dat = dat,dat_type = "local",R_matrix = R_matrix,
                       df_info = df_info,pval_threshold = pval_threshold,
                       F_threshold = F_threshold,
                       effect_size_cutoff = effect_size_cutoff,
                       extra_traits = extra_traits,
                       Filter = FALSE,
                       type_outcome = type_outcome, prevalence_outcome = prevalence_outcome,
                       prevalence_exposure = prevalence_exposure)
if (extra_traits == "None") {
  if (all(res$F.statistic > F_threshold)) {
    res_filter <- res
    stop_flag <- 'STOP'
  } else {
    row_to_remove <- which.min(res$F.statistic)
    res_filter <- res[-row_to_remove, ]
    stop_flag <- 'CONTINUE'
  }
} else {
  extra_trait_row <- which(res$id == extra_traits)
  if (all(res$F.statistic[-extra_trait_row] > F_threshold)) {
    res_filter <- res
    stop_flag <- 'STOP'
  } else {
    row_to_remove <- which.min(res$F.statistic[-extra_trait_row])
    if (row_to_remove >= extra_trait_row) {
      row_to_remove <- row_to_remove + 1
    }
    res_filter <- res[-row_to_remove, ]
    stop_flag <- 'CONTINUE'
  }
}
write(stop_flag, file = out_flag)
write.csv(data.frame(id = res_filter$id),file = out_id_list,row.names = F)
write.csv(res,file = out_df_strength,row.names = F)

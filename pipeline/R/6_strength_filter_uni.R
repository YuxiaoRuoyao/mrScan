library(dplyr)
library(mrScan)
library(MVMR)

dat_files <- unlist(snakemake@input[["dat"]])
R_files <- unlist(snakemake@input[["R"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
extra_traits <- snakemake@params[["extra_traits"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_df_strength <- snakemake@output[["out_df_strength"]]

all_res <- list()
for (i in seq_along(dat_files)) {
  dat <- readRDS(dat_files[i])
  R_matrix <- as.matrix(readRDS(R_files[i])$Re)
  res <- strength_filter(dat = dat, R_matrix = R_matrix, df_info = df_info,
                                  pval_threshold = pval_threshold,
                                  F_threshold = F_threshold,
                                  extra_traits = extra_traits,
                                  Filter = TRUE)
  res <- cbind(res[2,],res[1,])
  colnames(res)[4:6] <- c("F_main_exp","id_main_exp","trait_main_exp")
  all_res[[i]] <- res
}
final_res <- bind_rows(all_res)
filter_id <- all_res %>%
  filter(F.statistic > F_threshold & F_main_exp > F_threshold) %>% pull(id)

write.csv(data.frame(id = filter_id),file = out_id_list,row.names = F)
write.csv(final_res,file = out_df_strength,row.names = F)

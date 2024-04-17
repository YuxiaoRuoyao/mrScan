library(dplyr)
library(mrScan)
library(MVMR)

dat_files <- unlist(snakemake@input[["dat"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
id_exposure <- snakemake@params[["id_exposure"]]
extra_trait <- snakemake@params[["extra_trait"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]
out_df_strength <- snakemake@output[["out_df_strength"]]

all_res <- list()
for (i in seq_along(dat_files)) {
  dat <- readRDS(dat_files[i])
  res <- strength_filter(dat = dat, dat_type = "IEU",df_info = df_info,
                         pval_threshold = pval_threshold,
                         F_threshold = F_threshold,
                         Filter = FALSE)
  variable_id_row <- which(!(res$id %in% c(id_exposure, extra_trait)))
  res <- res[c(variable_id_row, which(res$id == id_exposure), which(res$id == extra_trait)), ]
  df <- data.frame(t(c(as.matrix(res))))
  colnames(df) <- paste0(rep(colnames(res), each = nrow(res)), 1:nrow(res))
  all_res[[i]] <- df
}
final_res <- bind_rows(all_res) %>%
  mutate(F.statistic1 = as.numeric(F.statistic1),
         F.statistic2 = as.numeric(F.statistic2),
         F.statistic3 = as.numeric(F.statistic3))
filter_id <- final_res %>%
  filter(F.statistic1 > F_threshold & F.statistic2 > F_threshold & F.statistic3 > F_threshold) %>%
  pull(id1)
all_id <- final_res$id1
other.id <- all_id[!all_id %in% filter_id]
df_info[df_info$id %in% other.id,"status"] <- "Delete due to weak instruments strength"

write.csv(data.frame(id = filter_id),file = out_id_list,row.names = F)
write.csv(final_res,file = out_df_strength,row.names = F)
write.csv(df_info,file = out_trait_info,row.names = F)

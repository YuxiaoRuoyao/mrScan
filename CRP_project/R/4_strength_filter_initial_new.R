library(dplyr)
library(MVMR)

strength_files <- unlist(snakemake@input[["strength_files"]])
df_info <- read.csv(snakemake@input[["trait_info"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
extra_trait <- snakemake@params[["extra_trait"]]
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]
out_df_strength <- snakemake@output[["out_df_strength"]]

all_res <- lapply(strength_files, readRDS)
final_res <- bind_rows(all_res)
if (extra_trait != "None") {
  final_res <- final_res %>%
    mutate(F.statistic1 = as.numeric(F.statistic1),
           F.statistic2 = as.numeric(F.statistic2),
           F.statistic3 = as.numeric(F.statistic3))
  
  filter_id <- final_res %>%
    filter(F.statistic1 > F_threshold & F.statistic2 > F_threshold & F.statistic3 > F_threshold) %>%
    pull(id1)
} else {
  final_res <- final_res %>%
    mutate(F.statistic1 = as.numeric(F.statistic1),
           F.statistic2 = as.numeric(F.statistic2))
  filter_id <- final_res %>%
    filter(F.statistic1 > F_threshold & F.statistic2 > F_threshold) %>%
    pull(id1)
}
all_id <- final_res$id1
other.id <- all_id[!all_id %in% filter_id]
df_info[df_info$id %in% other.id, "status"] <- "Delete due to weak instruments strength"
write.csv(data.frame(id = filter_id),file = out_id_list,row.names = F)
write.csv(final_res,file = out_df_strength,row.names = F)
write.csv(df_info,file = out_trait_info,row.names = F)

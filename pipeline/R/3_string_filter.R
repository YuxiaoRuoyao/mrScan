library(dplyr)
library(mrScan)

id.list <- read.csv(snakemake@input[["id_list"]])$id
df_info <- read.csv(snakemake@input[["trait_info"]])
inst_files <- unlist(snakemake@input[["files"]])
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
extra_traits <- unlist(snakemake@params[["extra_traits"]])
out_id_list <- snakemake@output[["out_id_list"]]
out_trait_info <- snakemake@output[["out_trait_info"]]

file_ids <- sub(".*inst_(.*)\\.RDS$", "\\1", basename(inst_files))
named_inst_files <- setNames(inst_files, file_ids)
get_row_count <- function(id) {
  data <- readRDS(named_inst_files[id])
  data.frame(id = id, n = nrow(data))
}
df_inst_counts <- purrr::map_dfr(id.list, get_row_count)
res_filter <- string_filter(id.list = id.list, df_info = df_info,
                            df_inst_counts = df_inst_counts,
                            R2_cutoff = R2_cutoff, extra_traits = extra_traits)

write.csv(data.frame(id = res_filter$id.list),file = out_id_list,row.names = F)
write.csv(res_filter$trait.info,file = out_trait_info,row.names = F)

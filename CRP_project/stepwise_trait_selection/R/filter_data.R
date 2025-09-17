library(dplyr)
library(stringr)

dat <- readRDS(snakemake@input[["beta"]])
id_list <- read.csv(snakemake@input[["id_list"]])$id
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

id_list <- unique(c(id_outcome, id_exposure, id_list))
beta_hat <- dat %>% select(ends_with(".beta"))
nms <- str_replace(names(beta_hat), ".beta", "")
out_id <- nms[!nms %in% id_list]
patterns_to_exclude <- paste0(out_id, "\\..*")
columns_to_exclude <- names(dat)[str_detect(names(dat), 
                                 pattern = str_c(patterns_to_exclude, collapse = "|"))]
dat_filtered <- dat %>% select(-all_of(columns_to_exclude)) %>% as_tibble()
saveRDS(dat_filtered, file = out)
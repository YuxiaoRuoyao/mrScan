library(dplyr)
library(purrr)
library(gridExtra)

res <- unlist(snakemake@input[["mvmr_file"]])
prefix <- snakemake@params[["prefix"]]
downstream_method <- snakemake@params[["method"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

prefix <- sub("_$", "", prefix)
extract_fdr <- function(filename) {
  match <- regexec("_FDR_p_([0-9.e-]+)_", filename)
  match_data <- regmatches(filename, match)
  if (length(match_data[[1]]) > 1) {
    return(match_data[[1]][2])
  } else {
    return(NA)
  }
}
if(grepl("MVMR",downstream_method)){
  type <- res %>% 
    strsplit(paste0("results_stepwise_plus/",prefix)) %>% sapply(tail, 1) %>%
    strsplit("_MVMR_") %>% sapply(head, 1) %>%
    gsub('selection_', '', .) %>% strsplit("_seed") %>%
    sapply(head, 1) %>% sub("^_", "", .) %>% unique()
}else{
  type <- res %>% 
    strsplit(paste0("results_stepwise_plus/",prefix)) %>% sapply(tail, 1) %>%
    strsplit("_MVMR_") %>% sapply(head, 1) %>% strsplit("_MR_") %>% sapply(head, 1) %>%
    gsub('selection_', '', .) %>% strsplit("_seed") %>%
    sapply(head, 1) %>% sub("^_", "", .) %>% unique()
}
MVMR_methods <- sub(".*_(MVMR_.*)\\.RDS$", "\\1", res) %>% unique()
fdr_values <- res %>% sapply(extract_fdr) %>% unique()
all_res <- data.frame()
for (m in MVMR_methods) {
  for (t in type) {
    if (t == "literature" | t == "marginal_p_0.05") {
      type_name <- paste0("selection_", t)
    } else {
      type_name <- t
    }
    for (fdr in fdr_values) {
      matching_files <- grep(paste0(type_name, "_",downstream_method, "_FDR_p_", fdr, ".*_", m, "\\.RDS$"), res, value = TRUE)
      for (filename in matching_files) {
        plus_id <- stringr::str_match(filename, "plus_([^_]+)_MVMR_")[,2]
        sub_res <- readRDS(filename)$res.summary
        sub_res$type <- t
        sub_res$FDR <- fdr
        sub_res$plus <- plus_id
        all_res <- bind_rows(all_res, sub_res)
      }
    }
  }
}
all_res <- all_res %>% mutate(CI_lower=b-qnorm(0.975)*se, CI_higher=b + qnorm(0.975)*se) %>%
           mutate(odds=exp(b),CI_lower=exp(CI_lower),CI_higher=exp(CI_higher))
all_res[all_res$type == "final","type"] <- "plus"
write.csv(all_res,file = out,row.names=FALSE)

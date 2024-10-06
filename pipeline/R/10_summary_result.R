library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)

res_mvmr <- unlist(snakemake@input[["mvmr_file"]])
res_uvmr <- unlist(snakemake@input[["uvmr_file"]])
prefix <- snakemake@params[["prefix"]]
downstream_method <- snakemake@params[["method"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out1 <- snakemake@output[["out1"]]
out2 <- snakemake@output[["out2"]]

prefix <- sub("_$", "", prefix)
res_mvmr <- Filter(function(f) !is.null(readRDS(f)),res_mvmr)
res <- c(res_mvmr,res_uvmr)
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
    strsplit(paste0("results/",prefix)) %>% sapply(tail, 1) %>%
    strsplit("_MVMR_") %>% sapply(head, 1) %>%
    gsub('selection_', '', .) %>% strsplit("_seed") %>%
    sapply(head, 1) %>% sub("^_", "", .) %>% unique()
}else{
  type <- res %>%
    strsplit(paste0("results/",prefix)) %>% sapply(tail, 1) %>%
    strsplit("_MVMR_") %>% sapply(head, 1) %>% strsplit("_MR_") %>% sapply(head, 1) %>%
    gsub('selection_', '', .) %>% strsplit("_seed") %>%
    sapply(head, 1) %>% sub("^_", "", .) %>% unique()
}
MVMR_methods <- sub(".*_(MVMR_.*)\\.RDS$", "\\1", res) %>% unique()
fdr_values <- res_mvmr %>% sapply(extract_fdr) %>% unique()
all_res <- data.frame()
for (m in MVMR_methods) {
  for (t in type) {
    if (t == "literature" | t == "marginal_p_0.05") {
      type_name <- paste0("selection_", t)
    } else {
      type_name <- t
    }
    for (fdr in fdr_values) {
      matching_files <- grep(paste0(type_name, "_",downstream_method, "_FDR_p_", fdr, ".*_", m, "\\.RDS$"), res_mvmr, value = TRUE)
      for (filename in matching_files) {
        sub_res <- readRDS(filename)$res.summary
        sub_res$type <- t
        sub_res$FDR <- fdr
        all_res <- bind_rows(all_res, sub_res)
      }
    }
  }
}
for (m in MVMR_methods) {
  for (fdr in fdr_values) {
    matching_files <-  grep(paste0("selection_UVMR", ".*_", m, "\\.RDS$"), res_uvmr, value = TRUE)
    for (filename in matching_files) {
      sub_res <- readRDS(filename)$res.summary
      sub_res$type <- "UVMR"
      all_res <- bind_rows(all_res, sub_res)
    }
  }
}
all_res <- all_res %>% mutate(CI_lower=b-qnorm(0.975)*se, CI_higher=b + qnorm(0.975)*se) %>%
  mutate(odds=exp(b),CI_lower=exp(CI_lower),CI_higher=exp(CI_higher))
#all_res[all_res$type == "unique_traits","type"] <- "All"
all_res[all_res$type == "final","type"] <- "Stepwise"
plt<- all_res %>% filter(exposure==id_exposure) %>%
  filter(converge == TRUE | is.na(converge)) %>%
  filter(se < 1) %>%
  ggplot() +
  geom_vline(xintercept = 1) +
  geom_point(aes(y = type, x = odds, color = method,  group = method),
             position=position_dodge(width = 0.9), size = 3) +
  geom_errorbar(aes(y = type, xmin =CI_lower, xmax = CI_higher, color = method),
                position=position_dodge(width = 0.9)) +
  xlab("Odds Ratio (95% CI)") + coord_flip() +
  theme_bw() + ggtitle(paste0("Direct causal effect of ",id_exposure,"  on ",id_outcome))+
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 10),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom")
ggsave(out2, plot = plt,width = 20, height = 10)
write.csv(all_res,file = out1,row.names=FALSE)

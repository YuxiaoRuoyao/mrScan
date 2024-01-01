library(dplyr)

res <- readRDS(snakemake@input[["file"]])
res_bidirection <- readRDS(snakemake@input[["file_bidirection"]])
p_cutoff <- as.numeric(snakemake@params[["p_cutoff"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
df_bidirection <- res_bidirection$df_bidirection

df_bidirection <- df_bidirection %>%
  mutate(p_ZtoX = 2*(1-pnorm(abs(b_ZtoX/se_ZtoX))),
         p_ZtoY = 2*(1-pnorm(abs(b_ZtoY/se_ZtoY))))
id.select <- df_bidirection %>% filter(id %in% id.list) %>%
  filter(p_ZtoX < p_cutoff & p_ZtoY < p_cutoff) %>% pull(id)
df_info[df_info$id %in% id.select,"status"] <- "Select by marginal selection"

saveRDS(list(id.list=id.select,trait.info=df_info),file = out)

library(VariantAnnotation)
library(gwasvcf)
library(rlang)
library(readr)
library(purrr)
library(stringr)
library(dplyr)
library(ieugwasr)
library(mrScan)
library(sumstatFactors)

#source("R/helpers.R")

res <- readRDS(snakemake@input[["file"]])
df_download <- read.csv(snakemake@input[["download"]],header=F)
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
out <- snakemake@output[["out"]]

id_list <- res$id.list
df_info <- res$trait.info

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "-GCST.*\\.h\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)

fulldat <- format_combine_gwas(df_file = file_list, c = c, df_info = df_info)

saveRDS(fulldat, file=out)

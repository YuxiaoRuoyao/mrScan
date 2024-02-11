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


res_name <- snakemake@input[["file"]]
df_download <- read.csv(snakemake@input[["download"]],header=F)
df_info_exposure_outcome <- read.csv(snakemake@input[["file_info_exposure_outcome"]])
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

if (file.size(res_name) != 0) {
    res <- readRDS(res_name)
    id_list <- res$id.list
    id_list <- c(id_outcome,id_exposure,id_list)
    df_info <- res$trait.info %>% full_join(df_info_exposure_outcome)
} else {
    id_list <- c(id_outcome,id_exposure)
    df_info <- df_info_exposure_outcome
}

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "-GCST.*\\.h\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)
df <- data.frame(id=id_list) %>% left_join(file_list)

fulldat <- format_combine_gwas(df_file = df, c = c, df_info = df_info)
saveRDS(fulldat, file=out)

library(VariantAnnotation)
library(gwasvcf)
library(rlang)
library(readr)
library(purrr)
library(stringr)
library(dplyr)
library(ieugwasr)
library(mrScan)
library(GFA)

input_files <- snakemake@input[["file"]]
df_download <- read.csv(snakemake@input[["download"]],header=F)
df_info_exposure_outcome <- read.csv(snakemake@input[["file_info_exposure_outcome"]])
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

if (length(input_files) == 3 && !is.na(input_files[[3]])) {
  id_list <- read.csv(input_files[[1]])$id
  id_list <- unique(c(id_outcome, id_exposure, id_list))
  df_info <- read.csv(input_files[[2]]) %>%
    full_join(df_info_exposure_outcome) %>%
    distinct(id, sample_size, .keep_all = TRUE)
  df_extra <- read.csv(input_files[[3]])
} else if (length(input_files) == 2 || (length(input_files) == 3 && is.na(input_files[[3]]))) {
  id_list <- read.csv(input_files[[1]])$id
  id_list <- unique(c(id_outcome, id_exposure, id_list))
  df_info <- read.csv(input_files[[2]]) %>%
    full_join(df_info_exposure_outcome) %>%
    distinct(id, sample_size, .keep_all = TRUE)
  df_extra <- NULL
} else {
  id_list <- c(id_outcome, id_exposure)
  df_info <- df_info_exposure_outcome
  df_extra <- NULL
}

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "-GCST.*\\.h\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)
if(!is.null(df_extra)){
  file_list <- full_join(file_list,df_extra) %>% distinct(id,.keep_all = TRUE)
}
df <- data.frame(id=id_list) %>% left_join(file_list)

fulldat <- format_combine_gwas(df_file = df, c = c, df_info = df_info)
saveRDS(fulldat, file=out)

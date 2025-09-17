library(dplyr)
library(mrScan)
library(stringr)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
trait_id <- snakemake@params[["trait_id"]]
trait_info <- read.csv(snakemake@input[["trait_info"]])
file_path <- snakemake@params[["path"]]
df_info_exposure_outcome <- read.csv(snakemake@input[["file_info_exposure_outcome"]])
df_download <- read.csv(snakemake@input[["download"]],header=F)
r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
clump_kb <- snakemake@params[["clump_kb"]]
ref_path  <- snakemake@params[["ref_path"]]
type <- snakemake@params[["ld_prioritization"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
ld_files <- unlist(snakemake@input[["l2"]])
m_files <- unlist(snakemake@input[["m"]])
out_R <- snakemake@output[["out_R"]]
out_dat <- snakemake@output[["out_dat"]]

id_list <- unique(c(id_outcome,id_exposure,trait_id))
df_info <- trait_info %>% full_join(df_info_exposure_outcome) %>%
  distinct(id, sample_size, .keep_all = TRUE)

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "-GCST.*\\.h\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)
df_file <- data.frame(id=id_list) %>% left_join(file_list)
res <- format_ldsc_prune(df_file = df_file,df_info = df_info,r2_thresh = r2_thresh,
                         clump_kb = clump_kb,type = type,pthresh = pthresh,ref_path = ref_path,
                         ld_files = ld_files,m_files = m_files)

saveRDS(res$R,file = out_R)
saveRDS(res$dat,file = out_dat)

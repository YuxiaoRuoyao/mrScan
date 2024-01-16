library(VariantAnnotation)
library(gwasvcf)
library(rlang)
library(readr)
library(purrr)
library(stringr)
library(dplyr)
library(ieugwasr)

source("R/helpers.R")

res_name <- snakemake@input[["file"]]
df_download <- read.csv(snakemake@input[["download"]],header=F)
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

if (file.size(res_name) != 0) {
    res <- readRDS(res_name)
    id_list <- res$id.list
    id_list <- c(id_outcome,id_exposure,id_list)
    df_info <- res$trait.info
} else {
    id_list <- c(id_outcome,id_exposure)
    #df_info <- gwasinfo(id_outcome,id_exposure) #need to edit it later
}

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "-GCST.*\\.h\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)
df <- data.frame(id=id_list) %>% left_join(file_list)

fulldat <- map(seq(nrow(df)), function(i){
  f <- df$location[i]
  if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
    dat <- format_ieu_chrom(f, c)
  }else if(str_ends(f, "h.tsv.gz")){
    dat <- format_flat_chrom(f, c,
                             snp_name = "hm_rsid",
                             pos_name = "hm_pos",
                             chrom_name = "hm_chrom",
                             A1_name = "hm_effect_allele",
                             A2_name = "hm_other_allele",
                             beta_hat_name = "hm_beta",
                             se_name = "standard_error",
                             p_value_name = "p_value",
                             af_name = "hm_effect_allele_frequency",
                             sample_size_name = NA,
                             effect_is_or = FALSE)
    if(all(is.na(dat$sample_size))){
      if(df$id[i] == "ebi-a-GCST90029070"){
        dat$sample_size <- 575531
      } # need to delete it later
      #if(df$id[i] %in% df_info$id){
      #  dat$sample_size <- df_info %>% filter(id == df$id[i]) %>% pull(sample_size)
      #}#else{
       # dat$sample_size <- gwasinfo(df$id[i])$sample_size
      #}
    }
  }
  n <- df$id[i]
  pos_name <- as_name(paste0(n, ".pos"))
  beta_name <- as_name(paste0(n, ".beta"))
  se_name <- as_name(paste0(n, ".se"))
  p_name <- as_name(paste0(n, ".p"))
  z_name <- as_name(paste0(n, ".z"))
  ss_name <- as_name(paste0(n, ".ss"))

  #dat$sample_size[is.na(dat$sample_size)] <- df_info[df_info$id == n,"sample_size"]
  dat <-dat %>% dplyr::mutate(Z = beta_hat/se) %>%
    dplyr::rename(REF = A2, ALT = A1) %>%
    dplyr::select(chrom, snp, REF, ALT,
                  !!pos_name := pos,
                  !!beta_name := beta_hat,
                  !!se_name := se,
                  !!p_name := p_value,
                  !!z_name := Z,
                  !!ss_name := sample_size)
}) %>%
  purrr::reduce(full_join, by = c("chrom", "snp", "REF", "ALT"))
dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
  fulldat <- filter(fulldat, !snp %in% dup_snps)
}
saveRDS(fulldat, file=out)

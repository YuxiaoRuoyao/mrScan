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

id_list <- read.csv(snakemake@input[["id_list"]])$id
trait_info <- read.csv(snakemake@input[["trait_info"]])
df_download <- read.csv(snakemake@input[["download"]],header=F)
df_info_exposure_outcome <- read.csv(snakemake@input[["file_info_exposure_outcome"]])
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
out <- snakemake@output[["out"]]

format_flat_chrom_edit <- function(file, chrom,
                                   snp_name,
                                   pos_name,
                                   chrom_name,
                                   A1_name, A2_name,
                                   beta_hat_name,
                                   se_name,
                                   p_value_name,
                                   af_name,
                                   sample_size_name,
                                   effect_is_or
){
  if(!p_value_name %in% c("", "NA", NA)){
    pstring <- paste0(", `", p_value_name, "`='d'")
  }else{
    pstring <- ""
    p_value_name <- NA
  }
  if(!sample_size_name %in% c("", "NA", NA)){
    sstring <- paste0(", `", sample_size_name, "`='d'")
  }else{
    sstring <- ""
    sample_size_name <- NA
  }
  if(!pos_name %in% c("", "NA", NA)){
    posstring <- paste0(", `", pos_name, "`='d'")
  }else{
    posstring <- ""
    pos_name <- NA
  }
  if(!af_name %in% c("", "NA", NA)){
    afstring <- paste0(", `", af_name, "`='d'")
  }else{
    afstring <- ""
    af_name <- NA
  }
  col_string <- paste0("cols_only(`", snp_name, "`='c', `",
                       A1_name , "`='c', `", A2_name, "`='c', `",
                       beta_hat_name , "`='d', `", se_name, "`='d', `",
                       chrom_name, "`='c' ", posstring,
                       pstring,  sstring, afstring, ")")
  if(str_ends(file, "gz") ){
    h <- read_tsv(pipe(paste0("gzip -cd ", file, " | head -2")))
    n <- which(names(h) == chrom_name)
    awk_cmd <- paste0("gzip -cd ", file, " | awk '{if ($", n, " == \"", chrom, "\") print $0}' - ")
  }else{
    h <- read_tsv(pipe(paste0("head -2 ", file)))
    n <- which(names(h) == chrom_name)
    awk_cmd <- paste0("awk '{if ($", n, " == \"", chrom, "\") print $0}' ", file)
  }
  X <- read_tsv(pipe(awk_cmd), col_types = eval(parse(text = col_string)), col_names = names(h))

  if(effect_is_or){
    X$beta <- log(X[[beta_hat_name]])
    beta_hat <- "beta"
  }
  dat <- GFA:::gwas_format(X, snp_name, beta_hat_name, se_name, A1_name,
                     A2_name, chrom_name, pos_name,
                     p_value = p_value_name,
                     allele_freq = af_name,
                     sample_size = sample_size_name,
                     compute_pval = TRUE)
  return(dat)
}

format_combine_gwas_edit <- function(df_file,c,df_info){
  fulldat <- purrr::map(seq(nrow(df_file)), function(i){
    f <- df_file$location[i]
    if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
      dat <- format_ieu_chrom(f, c)
    }else if(str_ends(f, ".h.tsv.gz")){
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
        if(df_file$id[i] %in% df_info$id){
          dat$sample_size <- df_info %>% filter(id == df_file$id[i]) %>% pull(sample_size)
        }else{
          dat$sample_size <- gwasinfo(df_file$id[i])$sample_size
        }
      }
    }else if (str_ends(f, ".tsv.gz") && !str_ends(f, ".h.tsv.gz")) {
      dat <- format_flat_chrom_edit(f, c,
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
        if(df_file$id[i] %in% df_info$id){
          dat$sample_size <- df_info %>% filter(id == df_file$id[i]) %>% pull(sample_size)
        }else{
          dat$sample_size <- gwasinfo(df_file$id[i])$sample_size
        }
      }
    }
    n <- df_file$id[i]
    pos_name <- as_name(paste0(n, ".pos"))
    af_name <- as_name(paste0(n, ".af"))
    beta_name <- as_name(paste0(n, ".beta"))
    se_name <- as_name(paste0(n, ".se"))
    p_name <- as_name(paste0(n, ".p"))
    z_name <- as_name(paste0(n, ".z"))
    ss_name <- as_name(paste0(n, ".ss"))

    dat$sample_size[is.na(dat$sample_size)] <- df_info[df_info$id == n,"sample_size"]
    dat <-dat %>% dplyr::mutate(Z = beta_hat/se) %>%
      dplyr::rename(REF = A2, ALT = A1) %>%
      dplyr::select(chrom, snp, REF, ALT,
                    !!pos_name := pos,
                    !!af_name := allele_freq,
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
  return(fulldat)
}

id_list <- unique(c(id_outcome,id_exposure,id_list))
df_info <- trait_info %>% full_join(df_info_exposure_outcome) %>%
    distinct(id, sample_size, .keep_all = TRUE)

file_list <- df_download$V1 %>% strsplit("/") %>% sapply(tail,1) %>% head(-1) %>%
  data.frame() %>% setNames("location") %>% filter(!str_detect(location, "\\.tbi$")) %>%
  mutate(id = case_when(
    str_detect(location, "\\.vcf\\.gz$") ~ gsub("\\.vcf\\.gz$", "", location),
    str_detect(location, "GCST.*\\.tsv\\.gz$") ~ paste0("ebi-a-", str_extract(location, "GCST[0-9]+"))
  )) %>% filter(id %in% id_list)
file_list$location <- paste0(file_path,file_list$location)
df <- data.frame(id=id_list) %>% left_join(file_list)

fulldat <- format_combine_gwas_edit(df_file = df, c = c, df_info = df_info)
saveRDS(fulldat, file=out)

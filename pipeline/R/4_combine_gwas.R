library(VariantAnnotation)
library(gwasvcf)
library(rlang)
library(readr)
library(purrr)
library(stringr)
library(dplyr)
library(ieugwasr)

source("R/helpers.R")

res <- readRDS(snakemake@input[["file"]])
c <- as.numeric(snakemake@wildcards[["chrom"]])
file_path <- snakemake@params[["path"]]
out <- snakemake@output[["out"]]

id_list <- res$id.list
df_info <- res$trait.info
file_list <- paste0(file_path,id_list,".vcf.gz")

fulldat <- map(seq(length(id_list)), function(i){
  f <- file_list[i]
  if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
    dat <- format_ieu_chrom(f, c)
  }else{
    dat <- format_flat_chrom(f, c,
                             snp[i],
                             pos[i],
                             chrom[i],
                             A1[i],
                             A2[i],
                             beta_hat[i],
                             se[i],
                             p_value[i],
                             af[i],
                             sample_size[i],
                             as.logical(effect_is_or[i]))
    if(all(is.na(dat$sample_size))){
      dat$sample_size <- pub_sample_size[i]
    }
  }
  n <- id_list[i]
  pos_name <- as_name(paste0(n, ".pos"))
  beta_name <- as_name(paste0(n, ".beta"))
  se_name <- as_name(paste0(n, ".se"))
  p_name <- as_name(paste0(n, ".p"))
  z_name <- as_name(paste0(n, ".z"))
  ss_name <- as_name(paste0(n, ".ss"))

  dat$sample_size[is.na(dat$sample_size)] <- df_info[df_info$id == id_list[i],"sample_size"]
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

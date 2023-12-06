#' @title Combine and harmonize local GWAS summary data with multiple traits
#' @param id_list GWAS ID list. The first should be the outcome ID.
#' @param file_list Directory path of formatted GWAS summary data. Default is in the current work directory.
#' @param out_dir Output data path. Default is in the current work directory.
#' @param prefix Name prefix for the output. Default = NULL
#' @returns Save one dataframe per chromosome with columns for SNP info
#'
#' @import stringr
#' @import readr
#' @import VariantAnnotation
#' @import gwasvcf
#' @import dplyr
#' @import rlang
#' @import purrr
#' @export

format_combine <- function(id_list,file_list,out_dir=NULL,prefix=NULL,
                           snp = NA, pos = NA, chrom = NA, A1 =  NA, A2 = NA,
                           beta_hat = NA, se = NA, p_value = NA, af = NA, sample_size = NA,
                           effect_is_or = "no"){
  for (c in seq(1:22)) {
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
      }
      n <- id_list[i]
      pos_name <- as_name(paste0(n, ".pos"))
      z_name <- as_name(paste0(n, ".z"))
      ss_name <- as_name(paste0(n, ".ss"))
      dat <-dat %>%  mutate(Z = beta_hat/se) %>%
        rename(REF = A2, ALT = A1) %>%
        select(chrom, snp, REF, ALT,
               !!pos_name := pos,
               !!z_name := Z,
               !!ss_name := sample_size)
    }) %>%
      purrr::reduce(full_join, by = c("chrom", "snp", "REF", "ALT"))
    dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
    if(length(dup_snps) > 0){
      fulldat <- filter(fulldat, !snp %in% dup_snps)
    }
    # Save table of how traits are missing each SNP for LD clumping
    miss <- fulldat %>%
      select(ends_with(".z")) %>%
      is.na(.) %>%
      rowSums(.)
    ix <- which(miss == 0)
    saveRDS(fulldat[ix,], file=paste0(out_dir,prefix,".zmat.",c,".RDS"))
  }
}

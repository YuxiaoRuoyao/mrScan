#' @title Combine and harmonize local GWAS summary data with multiple traits
#' @param df_file A dataframe contain "id" (trait id) and "location" (GWAS summary data location) columns.
#' @param c Numeric chromosome number
#' @param df_info Dataframe of trait info from previous steps containing sample size
#' @returns Save one dataframe per chromosome with columns for SNP info
#'
#' @import stringr
#' @import readr
#' @import gwasvcf
#' @import dplyr
#' @import rlang
#' @importFrom purrr reduce map
#' @rawNamespace import(VariantAnnotation, except = c(select,fixed))
#' @export
format_combine_gwas <- function(df_file,c,df_info){
  fulldat <- purrr::map(seq(nrow(df_file)), function(i){
    f <- df_file$location[i]
    if(str_ends(f, "vcf.gz") | str_ends(f, "vcf.bgz")){
      dat <- format_ieu_chrom(f, c)
    }else if(str_ends(f, ".h.tsv.gz")){
      header <- names(read_table(pipe(paste0("gzip -cd ", f, " | head -2"))))
      if (all(c("hm_rsid", "hm_pos", "hm_chrom", "hm_effect_allele",
                "hm_other_allele", "hm_beta", "standard_error",
                "p_value", "hm_effect_allele_frequency") %in% header)) {
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
      } else if (all(c("base_pair_location", "chromosome", "effect_allele",
                       "other_allele", "beta", "standard_error",
                       "p_value", "effect_allele_frequency") %in% header)) {
        dat <- format_flat_chrom(f, c,
                                 snp_name = "rsid",
                                 pos_name = "base_pair_location",
                                 chrom_name = "chromosome",
                                 A1_name = "effect_allele",
                                 A2_name = "other_allele",
                                 beta_hat_name = "beta",
                                 se_name = "standard_error",
                                 p_value_name = "p_value",
                                 af_name = "effect_allele_frequency",
                                 sample_size_name = NA,
                                 effect_is_or = FALSE)
      }else {
        stop(paste("Unrecognized column names in file:", f))
      }
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

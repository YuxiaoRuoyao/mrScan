#' @title Combine and harmonize GWAS summary data with multiple traits
#' @param id_list GWAS ID list. The first should be the outcome ID.
#' @param dir Directory path of formatted GWAS summary data. Default is in the current work directory.
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
combine_gwas_data <- function(id_list,dir=NULL,out_dir=NULL,prefix=NULL){
  for (chrom in seq(1:22)) {
    fulldat <- map(seq(length(id_list)),function(i){
      f <- paste0(dir, id_list[i], ".vcf.bgz")
      n <- id_list[i]
      v <- query_chrompos_file(paste0(chrom, ":1-536870911"), f)
      pos_name <- as_name(paste0(n, ".pos"))
      beta_name <- as_name(paste0(n, ".beta"))
      se_name <- as_name(paste0(n, ".se"))
      p_name <- as_name(paste0(n, ".p"))
      dat <- vcf_to_tibble(v) %>%
        dplyr::rename(snp=rsid) %>%
        dplyr::mutate(Z  = ES/SE,
                      P  = 2*pnorm(-abs(Z)))
      ii <- which(!is.na(dat$ES))
      dat <- dat %>%
        dplyr::select(chr:=seqnames,snp, REF, ALT,
                      !!pos_name := start,
                      !!beta_name := ES,
                      !!se_name := SE,
                      !!p_name := P)
      return(dat)
    }) %>%
      purrr::reduce(full_join, by = c("chr", "snp", "REF", "ALT"))
    dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
    if(length(dup_snps) > 0){
      fulldat <- filter(fulldat, !snp %in% dup_snps)
    }
    # Save table of how traits are missing each SNP for LD clumping
    miss <- fulldat %>%
      select(ends_with(".beta"))
    miss <- rowSums(is.na(miss))
    pmin <- fulldat %>%
      select(ends_with(".p"))
    pmin <- suppressWarnings(apply(pmin[,-1], 1, function(x){min(x, na.rm=TRUE)}))
    df <- data.frame(snp = fulldat$snp, pmin = pmin, miss = miss)
    ix <- which(miss <= 0)
    saveRDS(fulldat[ix,], file=paste0(out_dir,prefix,".beta.",chrom,".RDS"))
    saveRDS(df[ix,], file=paste0(out_dir,prefix,".info.",chrom,".RDS"))
  }
}



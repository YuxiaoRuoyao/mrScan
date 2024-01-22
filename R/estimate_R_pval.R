#' @title Estimate genetic correlation matrix by pvalue
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param p_thresh pvalue cutoff for filtering data
#' @returns A matrix with pairwise sample overlapping correlation between traits
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import sumstatFactors
#' @export
estimate_R_pval <- function(beta_files,p_thresh=0.05){
  X <- purrr::map_dfr(beta_files, readRDS)
  beta_hat <- X %>%
    select(ends_with(".beta")) %>%
    as.matrix()
  se_hat <- X %>%
    select(ends_with(".se")) %>%
    as.matrix()
  nms <- colnames(beta_hat) %>% stringr::str_replace(".beta$", "")
  Rpt <- R_pt(B_hat = beta_hat,
              S_hat = se_hat,
              p_val_thresh = p_thresh,
              return_cor = TRUE,
              make_well_conditioned = TRUE)
  colnames(Rpt) <- rownames(Rpt) <- nms
  return(Rpt)
}

#' @title MVMR analysis
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @param MVMR_method IVW, IVW_instrument_specific, MRBEE, GRAPPLE, ESMR
#' @param pleio_p_thresh P-value threshold for determining if a specific IV has
#' sufficient evidence of horizontal pleiotropy to be removed from causal estimation. Default=0
#' @param R Correlation matrix for the outcome and all exposures
#' @returns A dataframe with trait estimates for correponding methods
#'
#' @import TwoSampleMR
#' @import GRAPPLE
#' @import MRBEE
#' @import dplyr
#' @import purrr
#' @export
MVMR_analysis_local <- function(ld_prune_file_dir=NULL,prefix=NULL,R,MVMR_method,
                                p_thresh = 5e-8,pleio_p_thresh = 0){
  files <- paste0(ld_prune_file_dir,prefix,".beta.ldpruned.",seq(1,22),".RDS")
  X <- purrr::map_dfr(files, readRDS)
  beta_hat <- X %>% select(ends_with(".beta"))
  se <- X %>% select(ends_with(".se"))
  p <- X %>% select(ends_with(".p"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < p_thresh)
  i <- ncol(beta_hat)
  if(MVMR_method == "IVW"){
    hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                  exposure_pval = as.matrix(p[ix, 2:i]),
                  exposure_se = as.matrix(se[ix,2:i]),
                  outcome_beta = data.frame(beta_hat)[ix,1],
                  outcome_pval = data.frame(p)[ix,1],
                  outcome_se = data.frame(se)[ix,1],
                  expname = data.frame(id.exposure = nms, exposure = nms),
                  outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
    res <- TwoSampleMR::mv_multiple(hdat)$result
    res$method <- "IVW"
  }else if(MVMR_method == "IVW_instrument_specific"){
    hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                  exposure_pval = as.matrix(p[ix, 2:i]),
                  exposure_se = as.matrix(se[ix,2:i]),
                  outcome_beta = data.frame(beta_hat)[ix,1],
                  outcome_pval = data.frame(p)[ix,1],
                  outcome_se = data.frame(se)[ix,1],
                  expname = data.frame(id.exposure = nms, exposure = nms),
                  outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
    res <- TwoSampleMR::mv_multiple(hdat,instrument_specific = TRUE)$result
    res$method <- "IVW_instrument_specific"
  }else if(MVMR_method == "GRAPPLE"){
    res <- run_grapple(beta.exposure =  beta_hat[ix, 2:i],
                       beta.outcome = beta_hat[ix,1],
                       se.exposure = se[ix,2:i],
                       se.outcome = se[ix,1],
                       R = R)
  }else if(MVMR_method == "MRBEE"){
    res <- run_MRBEE(beta.exposure =  beta_hat[ix, 2:i],
                     beta.outcome = beta_hat[ix,1],
                     se.exposure = se[ix,2:i],
                     se.outcome = se[ix,1],
                     pleio_p_thresh = pleio_p_thresh, R = R)
  }
  return(res)
}

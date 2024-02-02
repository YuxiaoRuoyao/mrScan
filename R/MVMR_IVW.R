#' @title Use IVW to do MVMR analysis by locally data
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @returns A dataframe of result summary
#'
#' @import TwoSampleMR
#' @import dplyr
#' @importFrom purrr map_dfr
#' @export
MVMR_IVW <- function(beta_files,pval_threshold=5e-8){
  X <- purrr::map_dfr(beta_files, readRDS)
  beta_hat <- X %>% select(ends_with(".beta"))
  se <- X %>% select(ends_with(".se"))
  p <- X %>% select(ends_with(".p"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  i <- ncol(beta_hat)
  if(i>2){
    hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                  exposure_pval = as.matrix(p[ix, 2:i]),
                  exposure_se = as.matrix(se[ix,2:i]),
                  outcome_beta = data.frame(beta_hat)[ix,1],
                  outcome_pval = data.frame(p)[ix,1],
                  outcome_se = data.frame(se)[ix,1],
                  expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                  outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
    res_F <- mv_multiple(hdat)$result
    res_T <- mv_multiple(hdat,instrument_specific = TRUE)$result
    res_F$method <- paste0("IVW_",pval_threshold)
    res_T$method <- paste0("IVW_T_",pval_threshold)
    res <- rbind(res_F,res_T) %>% select(exposure,b,se,pval,method) %>%
      rename("pvalue" = "pval")
  }else{
    hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                  exposure_pval = as.matrix(p[ix, 2:i]),
                  exposure_se = as.matrix(se[ix,2:i]),
                  outcome_beta = data.frame(beta_hat)[ix,1],
                  outcome_pval = data.frame(p)[ix,1],
                  outcome_se = data.frame(se)[ix,1],
                  expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                  outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
    res_F <- mv_multiple(hdat)$result
    res_F$method <- paste0("IVW_",pval_threshold)
    res <- res_F %>% select(exposure,b,se,pval,method) %>%
      rename("pvalue" = "pval")
  }
  return(res)
}

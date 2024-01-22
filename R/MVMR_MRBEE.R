#' @title Use MRBEE to do MVMR analysis by locally data
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param pleio_p_thresh pvalue cutoff for determining if a specific IV has sufficient evidence of horizontal pleiotropy to be removed.
#' Default = 0
#' @returns A dataframe of result summary
#'
#' @import MRBEE
#' @import dplyr
#' @import purrr
#' @export
MVMR_MRBEE <- function(beta_files,R_matrix,pval_threshold = 5e-8,pleio_p_thresh = 0){
  X <- purrr::map_dfr(beta_files, readRDS)
  beta_hat <- X %>% select(ends_with(".beta"))
  se <- X %>% select(ends_with(".se"))
  p <- X %>% select(ends_with(".p"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  o <- match(colnames(R_matrix), nms)
  beta_hat <- data.frame(beta_hat[, o],check.names = F)
  se <- data.frame(se[, o],check.names = F)
  i <- ncol(beta_hat)
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  bT <- list(R = R_matrix, Ncor = Inf,
             EstHarm = beta_hat[ix,],
             SEHarm =  se[ix,])
  pD <- prepData(bT,verbose =FALSE)
  fit <- MRBEE.IMRP(pD, PleioPThreshold = pleio_p_thresh)
  res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                            b = fit$CausalEstimates[-1],
                            se = sqrt(diag(fit$VCovCausalEstimates))[-1])
  res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
  res.summary$method <- paste0("MRBEE_",pval_threshold,"_pleio_",pleio_p_thresh)
  return(res.summary)
}

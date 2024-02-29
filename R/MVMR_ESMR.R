#' @title Use ESMR to do MVMR analysis by locally data
#' @param dat A data frame of combined GWAS summary data after LD pruning. The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @returns A dataframe of result summary
#'
#' @import dplyr
#' @import esmr
#' @importFrom purrr map_dfr
#' @export
MVMR_ESMR <- function(dat,R_matrix,pval_threshold=5e-8){
  beta_hat <- dat %>% dplyr::select(ends_with(".beta"))
  se <- dat %>% dplyr::select(ends_with(".se"))
  p <- dat %>% dplyr::select(ends_with(".p"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  o <- match(colnames(R_matrix), nms)
  beta_hat <- data.frame(beta_hat[, o],check.names = F)
  se <- data.frame(se[, o],check.names = F)
  i <- ncol(beta_hat)
  fit <- esmr(beta_hat_Y = beta_hat[,1],
              se_Y = se[,1],
              beta_hat_X = beta_hat[,2:i],
              se_X = se[, 2:i],
              R = R_matrix,
              pval_thresh = pval_threshold)
  res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                            b = fit$beta$beta_m,
                            se = fit$beta$beta_s)
  res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
  res.summary$method <- paste0("ESMR_",pval_threshold)
  return(res.summary)
}

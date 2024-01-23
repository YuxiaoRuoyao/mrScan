#' @title Use GRAPPLE to do MVMR analysis by locally data
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 1e-5
#' @returns A dataframe of result summary
#'
#' @import GRAPPLE
#' @import dplyr
#' @import purrr
#' @import stringr
#' @export
MVMR_GRAPPLE <- function(beta_files,R_matrix,pval_threshold = 1e-5){
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
  grapple_dat <- data.frame(cbind(beta_hat, se))
  names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(i-1)),
                          "se_out", paste0("se_exp", 1:(i-1)))
  grapple_dat$selection_pvals <- apply(p[,-1, drop = F],1, min)
  res <- grappleRobustEst(data = grapple_dat,
                          plot.it =FALSE,
                          p.thres = pval_threshold,
                          cor.mat = R_matrix)
  res_warning <- tryCatch(grappleRobustEst(data = grapple_dat,
                          plot.it =FALSE,
                          p.thres = pval_threshold,
                          cor.mat = R_matrix),
                  warning = function(w) w)
  if(inherits(res_warning,"warning")){
    notConverge <- res_warning$message %>% str_detect("Did not converge")
  }else{
    notConverge <- FALSE
  }
  if(i > 2){
    res.summary <- data.frame(exposure=colnames(beta_hat)[-1],
                              b=res$beta.hat,
                              se=sqrt(diag(res$beta.var)),
                              pvalue=res$beta.p.value,
                              method = paste0("GRAPPLE_",pval_threshold),
                              converge = !notConverge)
  }else{
    res.summary <- data.frame(exposure=colnames(beta_hat)[-1],
                              b=res$beta.hat,
                              se=sqrt(res$beta.var),
                              pvalue=res$beta.p.value,
                              method = paste0("GRAPPLE_",pval_threshold),
                              converge = !notConverge)
  }
  return(res.summary)
}

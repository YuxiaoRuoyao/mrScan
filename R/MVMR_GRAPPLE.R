#' @title Use GRAPPLE to do MVMR analysis by locally data
#' @param dat A data frame of combined GWAS summary data after LD pruning. The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 1e-5
#' @returns A dataframe of result summary
#'
#' @import GRAPPLE
#' @import dplyr
#' @import stringr
#' @importFrom purrr map_dfr
#' @export
MVMR_GRAPPLE <- function(dat,R_matrix,pval_threshold = 1e-5){
  beta_hat <- dat %>% select(ends_with(".beta"))
  se <- dat %>% select(ends_with(".se"))
  p <- dat %>% select(ends_with(".p"))
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

#' @title Use MRBEE to do MVMR analysis by locally data
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @returns A dataframe of result summary
#'
#' @import MRBEE
#' @import dplyr
#' @importFrom purrr map_dfr
#' @export
MVMR_MRBEE <- function(beta_files,R_matrix,pval_threshold = 5e-8){
  X <- purrr::map_dfr(beta_files, readRDS)
  p <- X %>% select(ends_with(".p"))
  z <- X %>% select(ends_with(".z"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(p)<-names(z)<-nms
  o <- match(colnames(R_matrix), nms)
  z <- data.frame(z[, o],check.names = F)
  i <- ncol(beta_hat)
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  # Make the last one be outcome for R matrix
  R_matrix <- R_matrix[c(nms[-1],nms[1]),c(nms[-1],nms[1])]
  fit <- MRBEE.IMRP(by=z[ix,1],bX=as.matrix(z[ix,-1]),
                    byse=rep(1,length(ix)),
                    bXse=matrix(1,length(ix),i-1),
                    Rxy=R_matrix)
  res.summary <- data.frame(exposure = names(fit$theta),
                            b = fit$theta,
                            se = sqrt(diag(fit$covtheta)))
  res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
  res.summary$method <- paste0("MRBEE_",pval_threshold)
  rownames(res.summary) <- NULL
  return(res.summary)
}

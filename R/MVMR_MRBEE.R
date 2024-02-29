#' @title Use MRBEE to do MVMR analysis by locally data
#' @param dat A data frame of combined GWAS summary data after LD pruning. The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param pleio_threshold pvalue threshold in pleiotropy detection. Default = 0
#' @returns A dataframe of result summary
#'
#' @import MRBEE
#' @import dplyr
#' @importFrom purrr map_dfr
#' @export
MVMR_MRBEE <- function(dat,R_matrix,pval_threshold = 5e-8,pleio_threshold = 0){
  p <- dat %>% select(ends_with(".p"))
  z <- dat %>% select(ends_with(".z"))
  nms <- stringr::str_replace(names(z), ".z", "")
  names(p)<-names(z)<-nms
  o <- match(colnames(R_matrix), nms)
  z <- data.frame(z[, o],check.names = F)
  i <- ncol(z)
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  # Make the last one be outcome for R matrix
  R_matrix <- R_matrix[c(nms[-1],nms[1]),c(nms[-1],nms[1])]
  if(i>2){
    fit <- MRBEE.IMRP(by=z[ix,1],bX=as.matrix(z[ix,-1]),
                      byse=rep(1,length(ix)),
                      bXse=matrix(1,length(ix),i-1),
                      Rxy=R_matrix,
                      pv.thres = pleio_threshold, var.est = "variance")
    res.summary <- data.frame(exposure = names(fit$theta),
                              b = fit$theta,
                              se = sqrt(diag(fit$covtheta)))
  }else{
    fit <- MRBEE.IMRP.UV(by = z[ix,1],bx = z[ix,-1],
                         byse = rep(1,length(ix)),
                         bxse = rep(1,length(ix)),
                         Rxy=R_matrix,
                         pv.thres = pleio_threshold, var.est="variance")
    res.summary <- data.frame(exposure = nms[-1],
                              b = fit$theta,
                              se = sqrt(fit$vartheta))
  }
  res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
  res.summary$method <- paste0("MRBEE_",pval_threshold,"_pleio_",pleio_threshold)
  rownames(res.summary) <- NULL
  return(res.summary)
}

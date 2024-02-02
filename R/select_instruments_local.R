#' @title Select instruments from local GWAS summary data for MVMR analysis
#' @param beta_files Paths of combined local GWAS data after LD clumping
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @returns A list of harmonized data for the outcome (mvdat_y) and a list of harmonized data for the main exposure (mvdat_x)
#'
#' @import dplyr
#' @importFrom purrr map_dfr
#' @export
select_instruments_local <- function(beta_files,pval_threshold = 5e-8){
  X <- purrr::map_dfr(beta_files, readRDS)
  beta_hat <- X %>% select(ends_with(".beta"))
  se <- X %>% select(ends_with(".se"))
  p <- X %>% select(ends_with(".p"))
  p_x <- p[,-1]
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  pmin <- apply(p[,-1, drop = F], 1, min)
  pmin_x <- apply(p_x[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  ix_x <- which(pmin_x < pval_threshold)
  i <- ncol(beta_hat)
  mvdat_y <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                   exposure_pval = as.matrix(p[ix, 2:i]),
                   exposure_se = as.matrix(se[ix,2:i]),
                   outcome_beta = data.frame(beta_hat)[ix,1],
                   outcome_pval = data.frame(p)[ix,1],
                   outcome_se = data.frame(se)[ix,1],
                   expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                   outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
  mvdat_x <-  list(exposure_beta = as.matrix(beta_hat[ix_x, 3:i]),
                   exposure_pval = as.matrix(p[ix_x, 3:i]),
                   exposure_se = as.matrix(se[ix_x,3:i]),
                   outcome_beta = data.frame(beta_hat)[ix_x,2],
                   outcome_pval = data.frame(p)[ix_x,2],
                   outcome_se = data.frame(se)[ix_x,2],
                   expname = data.frame(id.exposure = nms[-c(1,2)], exposure = nms[-c(1,2)]),
                   outname = data.frame(id.outcome = nms[2], outcome = nms[2]))
  return(list(mvdat_x=mvdat_x,mvdat_y=mvdat_y))
}

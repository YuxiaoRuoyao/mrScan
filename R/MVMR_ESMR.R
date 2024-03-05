#' @title Use ESMR to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @returns A dataframe of result summary
#'
#' @import dplyr
#' @import esmr
#' @importFrom purrr map_dfc
#' @export
MVMR_ESMR <- function(dat,R_matrix,pval_threshold=5e-8,type){
  if(type == "local"){
    beta_hat <- dat %>% dplyr::select(ends_with(".beta"))
    se <- dat %>% dplyr::select(ends_with(".se"))
    p <- dat %>% dplyr::select(ends_with(".p"))
    nms <- stringr::str_replace(names(beta_hat), ".beta", "")
    names(beta_hat)<-names(se)<-names(p)<-nms
    o <- match(colnames(R_matrix), nms)
    beta_hat <- data.frame(beta_hat[, o],check.names = F)
    se <- data.frame(se[, o],check.names = F)
    i <- ncol(beta_hat)
    tryCatch({
      fit <- esmr(beta_hat_Y = beta_hat[,1],
                  se_Y = se[,1],
                  beta_hat_X = beta_hat[,2:i],
                  se_X = se[, 2:i],
                  R = R_matrix,
                  pval_thresh = pval_threshold)
      res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                                b = fit$beta$beta_m, se = fit$beta$beta_s) %>%
        mutate(pvalue = 2*pnorm(-abs(b/se)),
               method = paste0("ESMR_",pval_threshold))
    }, error = function(e){
      message("Error in MR_ESMR: ", e$message)
      res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                                b = NA, se = NA, pvalue = NA,
                                method = paste0("ESMR_",pval_threshold))
    })
  }
  if(type == "IEU"){
    id.exposure <- colnames(dat$exposure_beta)
    id.outcome <- dat$outname$id.outcome
    ss.exposure <- gwasinfo(id.exposure)$sample_size
    ss.outcome <- gwasinfo(id.outcome)$sample_size
    z.exposure<- dat$exposure_beta/dat$exposure_se
    z.norm.exposure <- sweep(z.exposure,2,sqrt(ss.exposure),`/`)
    names(ss.exposure) <- id.exposure
    se.norm.exposure <- map_dfc(ss.exposure, ~ rep(1/sqrt(.x),length.out = nrow(dat$exposure_se))) %>%
      as.matrix()
    tryCatch({
      fit <- esmr(beta_hat_Y = (dat$outcome_beta/dat$outcome_se)/sqrt(ss.outcome),
                  se_Y = rep(1/sqrt(ss.outcome),length(dat$outcome_se)),
                  beta_hat_X = z.norm.exposure,
                  se_X = se.norm.exposure)
      res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                                b = fit$beta$beta_m, se = fit$beta$beta_s) %>%
        mutate(pvalue = 2 * pnorm(-abs(b / se)), method = "MVMR_ESMR")
    }, error = function(e){
      message("Error in MVMR_ESMR: ", e$message)
      res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                                b = NA, se = NA, pvalue = NA,
                                method = "MVMR_ESMR")
    })
  }
  return(res.summary)
}

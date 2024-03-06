#' @title Use ESMR to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param ss.outcome A numeric number of sample size for the outcome. You can provide it when type = "IEU". Default = NULL
#' @returns A dataframe of result summary
#'
#' @import dplyr
#' @import esmr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_ESMR <- function(dat,R_matrix,pval_threshold=5e-8,type,
                      ss.exposure=NULL,ss.outcome=NULL){
  if(type == "local"){
    #beta_hat <- dat %>% dplyr::select(ends_with(".beta"))
    #se <- dat %>% dplyr::select(ends_with(".se"))
    z <- dat %>% select(ends_with(".z"))
    p <- dat %>% select(ends_with(".p"))
    ss <- dat %>% select(ends_with(".ss"))
    #nms <- stringr::str_replace(names(beta_hat), ".beta", "")
    nms <- stringr::str_replace(names(z), ".z", "")
    #names(beta_hat)<-names(se)<-names(p)<-nms
    names(z)<-names(p)<-names(ss)<-nms
    N <- apply(ss, 2, median, na.rm = TRUE)
    z.norm <- sweep(z,2,sqrt(N),`/`) %>% data.frame(check.names = F)
    se.norm <- purrr::map_dfc(ss, ~ rep(1/sqrt(.x),length.out = nrow(z))) %>%
      data.frame(check.names = F)
    o <- match(colnames(R_matrix), nms)
    #beta_hat <- data.frame(beta_hat[, o],check.names = F)
    #se <- data.frame(se[, o],check.names = F)
    z.norm <- z.norm[,o]
    se.norm <- se.norm[,o]
    #i <- ncol(beta_hat)
    i <- ncol(z.norm)
    res.summary <- data.frame(exposure = colnames(z.norm)[-1],
                              b = NA, se = NA, pvalue = NA,
                              method = paste0("ESMR_",pval_threshold))
    tryCatch({
      fit <- esmr(beta_hat_Y = z.norm[,1],
                  se_Y = se.norm[,1],
                  beta_hat_X = z.norm[,2:i],
                  se_X = se.norm[, 2:i],
                  R = R_matrix,
                  pval_thresh = pval_threshold)
      res.summary <- data.frame(exposure = colnames(z.norm)[-1],
                                b = fit$beta$beta_m, se = fit$beta$beta_s) %>%
        mutate(pvalue = 2*pnorm(-abs(b/se)),
               method = paste0("ESMR_",pval_threshold))
    }, error = function(e){
      message("Error in MVMR_ESMR: ", e$message)
    })
  }
  if(type == "IEU"){
    id.exposure <- colnames(dat$exposure_beta)
    id.outcome <- dat$outname$id.outcome
    if(is.null(ss.exposure)){
      ss.exposure <- gwasinfo(id.exposure)$sample_size
    }
    if(is.null(ss.outcome)){
      ss.outcome <- gwasinfo(id.outcome)$sample_size
    }
    z.exposure<- dat$exposure_beta/dat$exposure_se
    z.norm.exposure <- sweep(z.exposure,2,sqrt(ss.exposure),`/`)
    names(ss.exposure) <- id.exposure
    se.norm.exposure <- map_dfc(ss.exposure, ~ rep(1/sqrt(.x),length.out = nrow(dat$exposure_se))) %>%
      as.matrix()
    res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                              b = NA, se = NA, pvalue = NA,
                              method = "MVMR_ESMR")
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
    })
  }
  return(res.summary)
}

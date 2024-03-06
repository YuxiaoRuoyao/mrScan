#' @title Use IVW to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param ss.outcome A numeric number of sample size for the outcome. You can provide it when type = "IEU". Default = NULL
#' @returns A dataframe of result summary
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_IVW <- function(dat,pval_threshold=5e-8,type,
                     ss.exposure = NULL, ss.outcome = NULL){
  if(type == "local"){
    #beta_hat <- dat %>% select(ends_with(".beta"))
    #se <- dat %>% select(ends_with(".se"))
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
    pmin <- apply(p[,-1, drop = F], 1, min)
    ix <- which(pmin < pval_threshold)
    #i <- ncol(beta_hat)
    i <- ncol(z.norm)
    if(i>2){
      hdat <-  list(exposure_beta = as.matrix(z.norm[ix, 2:i]),
                    exposure_pval = as.matrix(p[ix, 2:i]),
                    exposure_se = as.matrix(se.norm[ix,2:i]),
                    outcome_beta = data.frame(z.norm)[ix,1],
                    outcome_pval = data.frame(p)[ix,1],
                    outcome_se = data.frame(se.norm)[ix,1],
                    expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                    outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
      res_F <- mv_multiple(hdat)$result
      res_T <- mv_multiple(hdat,instrument_specific = TRUE)$result
      res_F$method <- paste0("IVW_",pval_threshold)
      res_T$method <- paste0("IVW_T_",pval_threshold)
      res <- rbind(res_F,res_T) %>% select(exposure,b,se,pval,method) %>%
        rename("pvalue" = "pval")
    }else{
      hdat <-  list(exposure_beta = as.matrix(z.norm[ix, 2:i]),
                    exposure_pval = as.matrix(p[ix, 2:i]),
                    exposure_se = as.matrix(se.norm[ix,2:i]),
                    outcome_beta = data.frame(z.norm)[ix,1],
                    outcome_pval = data.frame(p)[ix,1],
                    outcome_se = data.frame(se.norm)[ix,1],
                    expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                    outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
      res_F <- mv_multiple(hdat)$result
      res_F$method <- paste0("IVW_",pval_threshold)
      res <- res_F %>% select(exposure,b,se,pval,method) %>%
        rename("pvalue" = "pval")
    }
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
    dat_norm <- list(exposure_beta = z.norm.exposure, exposure_pval = dat$exposure_pval,
                     exposure_se = se.norm.exposure,
                     outcome_beta = (dat$outcome_beta/dat$outcome_se)/sqrt(ss.outcome),
                     outcome_pval = dat$outcome_pval,
                     outcome_se = rep(1/sqrt(ss.outcome),length(dat$outcome_se)),
                     expname = dat$expname, outname = dat$outname)
    res <- mv_multiple(dat_norm)$result %>% select(id.exposure,id.outcome,b,se,pval) %>%
      rename("pvalue" = "pval") %>% mutate(method = "MVMR_IVW")
  }
  return(res)
}

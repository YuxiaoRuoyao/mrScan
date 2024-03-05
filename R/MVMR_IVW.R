#' @title Use IVW to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @returns A dataframe of result summary
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_IVW <- function(dat,pval_threshold=5e-8,type){
  if(type == "local"){
    beta_hat <- dat %>% select(ends_with(".beta"))
    se <- dat %>% select(ends_with(".se"))
    p <- dat %>% select(ends_with(".p"))
    nms <- stringr::str_replace(names(beta_hat), ".beta", "")
    names(beta_hat)<-names(se)<-names(p)<-nms
    pmin <- apply(p[,-1, drop = F], 1, min)
    ix <- which(pmin < pval_threshold)
    i <- ncol(beta_hat)
    if(i>2){
      hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                    exposure_pval = as.matrix(p[ix, 2:i]),
                    exposure_se = as.matrix(se[ix,2:i]),
                    outcome_beta = data.frame(beta_hat)[ix,1],
                    outcome_pval = data.frame(p)[ix,1],
                    outcome_se = data.frame(se)[ix,1],
                    expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                    outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
      res_F <- mv_multiple(hdat)$result
      res_T <- mv_multiple(hdat,instrument_specific = TRUE)$result
      res_F$method <- paste0("IVW_",pval_threshold)
      res_T$method <- paste0("IVW_T_",pval_threshold)
      res <- rbind(res_F,res_T) %>% select(exposure,b,se,pval,method) %>%
        rename("pvalue" = "pval")
    }else{
      hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
                    exposure_pval = as.matrix(p[ix, 2:i]),
                    exposure_se = as.matrix(se[ix,2:i]),
                    outcome_beta = data.frame(beta_hat)[ix,1],
                    outcome_pval = data.frame(p)[ix,1],
                    outcome_se = data.frame(se)[ix,1],
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
    ss.exposure <- gwasinfo(id.exposure)$sample_size
    ss.outcome <- gwasinfo(id.outcome)$sample_size
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

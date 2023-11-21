#' @title MVMR analysis
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param r2 LD-clump r2. Default = 0.001
#' @param kb LD-clump kb. Default = 10000
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @param find_proxies Whether look for proxies. Default = TRUE
#' @param harmonise_strictness Data harmonise strictness. See documentation of TwoSample MR
#' @param pop Super population to use. Default = "EUR"
#' @param MVMR_method IVW, IVW_instrument_specific, MVMR_median, MRBEE, GRAPPLE, ESMR
#' @param pleio_p_thresh P-value threshold for determining if a specific IV has
#' sufficient evidence of horizontal pleiotropy to be removed from causal estimation. Default=0
#' @param Rcor Correlation matrix for the outcome and all exposures
#' @returns A dataframe with trait estimates for correponding methods
#'
#' @import TwoSampleMR
#' @import GRAPPLE
#' @import MRBEE
#' @import MendelianRandomization
#' @export
MVMR_analysis <- function(id_exposure,id_outcome,id.list,df_info,r2 = 0.001, kb = 10000,
                          pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                          harmonise_strictness = 2, MVMR_method,
                          pleio_p_thresh = 0, Rcor = NULL){
  inst <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,
                                              find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out <- TwoSampleMR::extract_outcome_data(inst$SNP, id_outcome)
  mvdat <- TwoSampleMR::mv_harmonise_data(inst, out)
  if(MVMR_method == "IVW"){
    res <- TwoSampleMR::mv_multiple(mvdat)$result
    res$method <- "IVW"
  }else if(MVMR_method == "IVW_instrument_specific"){
    res <- TwoSampleMR::mv_multiple(mvdat,instrument_specific = TRUE)$result
    res$method <- "IVW_instrument_specific"
  }else if(MVMR_method == "MVMR_median"){
    mvmr_obj <- MendelianRandomization::mr_mvinput(bx = mvdat$exposure_beta,
                                       bxse = mvdat$exposure_se,
                                       by = mvdat$outcome_beta,
                                       byse = mvdat$outcome_se)
    fit <- MendelianRandomization::mr_mvmedian(mvmr_obj)
    res <- data.frame(exposure = colnames(mvdat$exposure_beta),
                      b = fit@Estimate,
                      se = fit@StdError,
                      pvalue = fit@Pvalue, method = "MVMR_median")

  }else if(MVMR_method == "GRAPPLE"){
    res <- run_grapple(beta.exposure =  mvdat$exposure_beta,
                       beta.outcome = mvdat$outcome_beta,
                       se.exposure = mvdat$exposure_se,
                       se.outcome = mvdat$outcome_se)
  }else if(MVMR_method == "MRBEE"){
    res <- run_MRBEE(beta.exposure =  mvdat$exposure_beta,
                     beta.outcome = mvdat$outcome_beta,
                     se.exposure = mvdat$exposure_se,
                     se.outcome = mvdat$outcome_se,
                     pleio_p_thresh = pleio_p_thresh, Rcor = Rcor)
  }
  # still need the ESMR function
  return(res)
}

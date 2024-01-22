#' @title Select instruments by api for MVMR analysis
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param r2 LD-clump r2. Default = 0.001
#' @param kb LD-clump kb. Default = 10000
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @param find_proxies Whether look for proxies. Default = TRUE
#' @param pop Super population to use. Default = "EUR"
#' @param harmonise_strictness Data harmonise strictness. See documentation of TwoSample MR
#' @returns A list of harmonized data for the outcome (mvdat_y) and a list of harmonized data for the main exposure (mvdat_x)
#'
#' @import TwoSampleMR
#' @export
select_instruments_api <- function(id.list,id_exposure,id_outcome,r2 = 0.001, kb = 10000,
                                   pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                                   harmonise_strictness = 2){
  # X ~ Z
  inst_x <- mv_extract_exposures(id.list,clump_r2 = r2, clump_kb = kb,
                                 harmonise_strictness = harmonise_strictness,find_proxies = find_proxies,
                                 pval_threshold = pval_threshold, pop = pop)
  out_x <- extract_outcome_data(inst_x$SNP, id_exposure)
  mvdat_x <- mv_harmonise_data(inst_x, out_x)
  # Y ~ Z
  inst_y <- mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                 harmonise_strictness = harmonise_strictness,
                                 find_proxies = find_proxies,
                                 pval_threshold = pval_threshold, pop = pop)
  out_y <- extract_outcome_data(inst_y$SNP, id_outcome)
  mvdat_y <- mv_harmonise_data(inst_y, out_y)
  return(list(mvdat_x=mvdat_x,mvdat_y=mvdat_y))
}

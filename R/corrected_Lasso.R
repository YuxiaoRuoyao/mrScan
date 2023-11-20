#' @title Use stepwise selection to do confounder selection
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
#' @param radius_type radius type for corresponding loss. Default = "1se"
#' @param seed Default = 1
#' @param maxits maximum number of iterations for algorithm converges. Default=1000000
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import TwoSampleMR
#' @import hdme
#' @import dplyr
#' @export
corrected_Lasso <- function(id_exposure,id_outcome,id.list,df_info,r2 = 0.001, kb = 10000,
                            pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                            harmonise_strictness = 2,
                            radius_type="1se",seed = 1,maxits = 1000000){
  if(radius_type == "min"){
    radius <- 'radius_min'
  }else{
    radius <- 'radius_1se'
  }
  inst_y <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,
                                              find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out_y <- TwoSampleMR::extract_outcome_data(inst_y$SNP, id_outcome)
  mvdat_y <- mv_harmonise_data(inst_y, out_y)
  ids <- colnames(mvdat_y$exposure_beta)
  ids <- ids[!ids %in% id_exposure]
  ss <- df_info %>% filter(id %in% ids) %>%
    arrange(match(id, ids)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_y <- hdme::cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                             y = mvdat_y$outcome_beta,
                                             sigmaUU = diag(1/ss))
  corrected_y <- hdme::corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                       y = mvdat_y$outcome_beta,
                                       sigmaUU = diag(1/ss),
                                       radii = cv_corrected_y[[radius]],
                                       maxits = maxits)
  coef_y <-coef(corrected_y)
  id.select <- ids[coef_y$coefficient]
  df_info[df_info$id %in% id.select,"status"] <- "Select by Corrected Lasso"
  return(list(id.list=id.select,trait.info=df_info))
}

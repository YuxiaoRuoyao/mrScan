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
double_corrected_Lasso <- function(id_exposure,id_outcome,id.list,df_info,r2 = 0.001, kb = 10000,
                            pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                            harmonise_strictness = 2,
                            radius_type="1se",seed = 1, maxits = 1000000){
  if(radius_type == "min"){
    radius <- 'radius_min'
  }else{
    radius <- 'radius_1se'
  }
  # X ~ Z
  inst_x <- TwoSampleMR::mv_extract_exposures(id.list,clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out_x <- TwoSampleMR::extract_outcome_data(inst_x$SNP, id_exposure)
  mvdat_x <- mv_harmonise_data(res_inst, out_x)
  ids_x <- colnames(mvdat_x$exposure_beta)
  ss_x <- df_info %>% filter(id %in% ids_x) %>%
    arrange(match(id, ids_x)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_x <- hdme::cv_corrected_lasso(W = mvdat_x$exposure_beta,
                                             y = mvdat_x$outcome_beta,
                                             sigmaUU = diag(1/ss_x))
  corrected_x <- hdme::corrected_lasso(W = mvdat_x$exposure_beta,
                                       y = mvdat_x$outcome_beta,
                                       sigmaUU = diag(1/ss_x),
                                       radii = cv_corrected_x[[radius]],
                                       maxits = maxits)
  coef_x <-coef(corrected_x)
  id_x <- ids_x[coef_x$coefficient]
  # Y ~ Z
  inst_y <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,
                                              find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out_y <- TwoSampleMR::extract_outcome_data(inst_y$SNP, id_outcome)
  mvdat_y <- mv_harmonise_data(inst_y, out_y)
  ids_y <- colnames(mvdat_y$exposure_beta)
  ids_y <- ids_y[!ids_y %in% id_exposure]
  ss_y <- df_info %>% filter(id %in% ids_y) %>%
    arrange(match(id, ids_y)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_y <- hdme::cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                             y = mvdat_y$outcome_beta,
                                             sigmaUU = diag(1/ss_y))
  corrected_y <- hdme::corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                       y = mvdat_y$outcome_beta,
                                       sigmaUU = diag(1/ss_y),
                                       radii = cv_corrected_y[[radius]],
                                       maxits = maxits)
  coef_y <-coef(corrected_y)
  id_y <- ids_y[coef_y$coefficient]
  id.select <- union(id_x,id_y)
  df_info[df_info$id %in% id.select,"status"] <- "Select by Double Corrected Lasso"
  return(list(id.list=id.select,trait.info=df_info))
}

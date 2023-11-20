#' @title Use classic Lasso to do confounder selection
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
#' @param lambda_type Lasso penalty type. Default = "min"
#' @param seed Default = 1
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import TwoSampleMR
#' @import glmnet
#' @export
classic_Lasso <- function(id_exposure, id_outcome, id.list,df_info,r2 = 0.001, kb = 10000,
                         pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                         harmonise_strictness = 2, lambda_type = "min",seed = 1){
  if(lambda_type == "min"){
    penalty <- 'lambda.min'
  }else{
    penalty <- 'lambda.1se'
  }
  # Y ~ Z
  inst_y <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,
                                              find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out_y <- TwoSampleMR::extract_outcome_data(inst_y$SNP, id_outcome)
  mvdat_y <- mv_harmonise_data(inst_y, out_y)
  set.seed(seed)
  # not penalize X
  cv_model_y <- glmnet::cv.glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                  y=mvdat_y$outcome_beta, alpha = 1)
  best_model_y <- glmnet::glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                 y=mvdat_y$outcome_beta, alpha = 1,
                                 lambda = cv_model_y[[penalty]])
  A_y<-rownames(coef(best_model_y))[which(coef(best_model_y)[,1]!=0)][-1]
  df_info[df_info$id %in% A_y,"status"] <- "Select by Classic Lasso"
  return(list(id.list=A_y,trait.info=df_info))
}

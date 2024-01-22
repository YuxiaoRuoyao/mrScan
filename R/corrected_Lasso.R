#' @title Use stepwise selection to do confounder selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param mvdat_y Harmonized data for the outcome produced by select_instruments step
#' @param radius_type radius type for corresponding loss. Default = "1se"
#' @param seed Default = 1
#' @param maxits maximum number of iterations for algorithm converges. Default=1000000
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import hdme
#' @import dplyr
#' @export
corrected_Lasso <- function(id_exposure,id.list,df_info,mvdat_y,
                            radius_type="1se",seed = 1,maxits = 1000000){
  if(radius_type == "min"){
    radius <- 'radius_min'
  }else{
    radius <- 'radius_1se'
  }
  ids <- colnames(mvdat_y$exposure_beta)
  ids <- ids[!ids %in% id_exposure]
  ss <- df_info %>% filter(id %in% ids) %>%
    arrange(match(id, ids)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_y <- cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                       y = mvdat_y$outcome_beta,
                                       sigmaUU = diag(1/ss))
  corrected_y <- corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                 y = mvdat_y$outcome_beta,
                                 sigmaUU = diag(1/ss),
                                 radii = cv_corrected_y[[radius]],
                                 maxits = maxits)
  coef_y <-coef(corrected_y)
  id.select <- ids[coef_y$coefficient]
  df_info[df_info$id %in% id.select,"status"] <- "Select by Corrected Lasso"
  return(list(id.list=id.select,trait.info=df_info))
}

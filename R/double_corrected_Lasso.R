#' @title Use stepwise selection to do confounder selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param mvdat_x Harmonized data for the main exposure produced by select_instruments step
#' @param mvdat_y Harmonized data for the outcome produced by select_instruments step
#' @param radius_type radius type for corresponding loss. Default = "1se"
#' @param seed Default = 1
#' @param maxits maximum number of iterations for algorithm converges. Default=1000000
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import hdme
#' @import dplyr
#' @export
double_corrected_Lasso <- function(id_exposure,id.list,df_info,
                                   mvdat_x,mvdat_y,
                                   radius_type="1se",seed = 1, maxits = 1000000){
  if(radius_type == "min"){
    radius <- 'radius_min'
  }else{
    radius <- 'radius_1se'
  }
  ids_x <- colnames(mvdat_x$exposure_beta)
  ss_x <- df_info %>% filter(id %in% ids_x) %>%
    arrange(match(id, ids_x)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_x <- cv_corrected_lasso(W = mvdat_x$exposure_beta,
                                       y = mvdat_x$outcome_beta,
                                       sigmaUU = diag(1/ss_x))
  corrected_x <- corrected_lasso(W = mvdat_x$exposure_beta,
                                 y = mvdat_x$outcome_beta,
                                 sigmaUU = diag(1/ss_x),
                                 radii = cv_corrected_x[[radius]],
                                 maxits = maxits)
  coef_x <-coef(corrected_x)
  id_x <- ids_x[coef_x$coefficient]

  ids_y <- colnames(mvdat_y$exposure_beta)
  ids_y <- ids_y[!ids_y %in% id_exposure]
  ss_y <- df_info %>% filter(id %in% ids_y) %>%
    arrange(match(id, ids_y)) %>% pull(sample_size)
  set.seed(seed)
  cv_corrected_y <- cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                       y = mvdat_y$outcome_beta,
                                       sigmaUU = diag(1/ss_y))
  corrected_y <- corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
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

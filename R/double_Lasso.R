#' @title Use classic Double Lasso to do confounder selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param mvdat_x Harmonized data for the main exposure produced by select_instruments step
#' @param mvdat_y Harmonized data for the outcome produced by select_instruments step
#' @param lambda_type Lasso penalty type. Default = "min"
#' @param seed Default = 1
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import TwoSampleMR
#' @import glmnet
#' @export
double_Lasso <- function(id_exposure, id.list,df_info,
                         mvdat_x,mvdat_y,
                         lambda_type = "min",seed = 1){
  if(lambda_type == "min"){
    penalty <- 'lambda.min'
  }else{
    penalty <- 'lambda.1se'
  }
  set.seed(seed)
  cv_model_x <- cv.glmnet(x=mvdat_x$exposure_beta, y=mvdat_x$outcome_beta,
                          alpha = 1)
  best_model_x <- glmnet(x=mvdat_x$exposure_beta,
                         y=mvdat_x$outcome_beta, alpha = 1,
                         lambda = cv_model_x[[penalty]])
  A_x<-rownames(coef(best_model_x))[which(coef(best_model_x)[,1]!=0)][-1]

  set.seed(seed)
  cv_model_y <- cv.glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                          y=mvdat_y$outcome_beta, alpha = 1)
  best_model_y <- glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                         y=mvdat_y$outcome_beta, alpha = 1,
                         lambda = cv_model_y[[penalty]])
  A_y<-rownames(coef(best_model_y))[which(coef(best_model_y)[,1]!=0)][-1]
  A<-union(A_x,A_y)
  df_info[df_info$id %in% A,"status"] <- "Select by Double Lasso"
  return(list(id.list=A,trait.info=df_info))
}

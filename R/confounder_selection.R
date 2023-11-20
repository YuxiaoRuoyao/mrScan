#' @title Select confounders by variable selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on downstream_filter
#' @param df_info Dataframe of trait info from previous steps
#' @param method Variable selection method. Include Lasso-type, marginal selection and self-input
#' @param stepwise_method either forward, backward or both. Default = "forward"
#' @param r2 LD-clump r2. Default = 0.001
#' @param kb LD-clump kb. Default = 10000
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @param find_proxies Whether look for proxies. Default = TRUE
#' @param pop Super population to use. Default = "EUR"
#' @param harmonise_strictness Data harmonise strictness. See documentation of TwoSample MR
#' @param lambda_type Lasso penalty type. Default = "min"
#' @param radius_type Radius type for corresponding loss to corrected Lasso. Default = "min"
#' @param seed Default = 1
#' @param maxits maximum number of iterations for algorithm converges. Default=1000000
#' @param p_cutoff pvalue threshold for marginal selection
#' @param df_bidirection dataframe of traits with estimates of four direction, used
#' by marginal selection. Results from downstream_filter
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import dplyr
#' @import data.table
#' @export
confounder_selection <- function(id_exposure,id_outcome,id.list,df_info,method,
                                 stepwise_method = "forward",
                                 r2 = 0.001, kb = 10000,
                                 pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                                 harmonise_strictness = 2, lambda_type = "min",
                                 radius_type = "1se",seed = 1,maxits = 1000000,
                                 p_cutoff = 0.05, df_bidirection = NULL){
  if(method == "classic_Lasso"){
    res <- classic_Lasso(id_exposure = id_exposure, id_outcome = id_outcome,
                                       id.list=id.list, df_info = df_info, r2=r2,
                                       kb=kb, pval_threshold=pval_threshold,
                                       find_proxies=find_proxies, pop=pop,
                                       harmonise_strictness=harmonise_strictness,
                                       lambda_type=lambda_type,seed=seed)
  }else if(method == "double_Lasso"){
    res <- double_Lasso(id_exposure = id_exposure, id_outcome = id_outcome,
                                     id.list=id.list, df_info = df_info, r2=r2,
                                     kb=kb, pval_threshold=pval_threshold,
                                     find_proxies=find_proxies, pop=pop,
                                     harmonise_strictness=harmonise_strictness,
                                     lambda_type=lambda_type,seed=seed)
  }else if(method == "stepwise"){
    res <- stepwise(id_exposure = id_exposure, id_outcome = id_outcome,
                             id.list = id.list, df_info = df_info,
                             method = stepwise_method, r2=r2, kb=kb,
                             pval_threshold=pval_threshold,find_proxies=find_proxies,
                             pop=pop,harmonise_strictness=harmonise_strictness)
  }else if(method == "marginal"){
    if(is.null(df_bidirection)){
      stop("Need to input dataframe of bidirection estimates!")
    }else{
      res <- marginal(id.list = id.list,df_info = df_info,
                               df_bidirection = df_bidirection,
                               p_cutoff = p_cutoff)
    }
  }else if(method == "corrected_Lasso"){
    res <- corrected_Lasso(id_exposure = id_exposure, id_outcome = id_outcome,
                                     id.list = id.list, df_info = df_info,r2=r2,
                                     kb=kb, pval_threshold=pval_threshold,
                                     find_proxies=find_proxies, pop=pop,
                                     harmonise_strictness=harmonise_strictness,
                                     radius_type="1se",seed=seed,maxits = maxits)
  }else if(method == "double_corrected_Lasso"){
    res <- double_corrected_Lasso(id_exposure = id_exposure, id_outcome = id_outcome,
                                                   id.list = id.list, df_info = df_info,r2=r2,
                                                   kb=kb, pval_threshold=pval_threshold,
                                                   find_proxies=find_proxies, pop=pop,
                                                   harmonise_strictness=harmonise_strictness,
                                                   radius_type="1se",seed=seed,maxits = maxits)
  }
  return(res)
}

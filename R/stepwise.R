#' @title Use stepwise selection to do confounder selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param method either forward, backward or both. Default = forward
#' @param r2 LD-clump r2. Default = 0.001
#' @param kb LD-clump kb. Default = 10000
#' @param pval_threshold Instrument selection cutoff. Default = 5e-8
#' @param find_proxies Whether look for proxies. Default = TRUE
#' @param harmonise_strictness Data harmonise strictness. See documentation of TwoSample MR
#' @param pop Super population to use. Default = "EUR"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import TwoSampleMR
#' @import stringr
#' @export
stepwise <- function(id_exposure, id_outcome, id.list,df_info,method = "forward",
                     r2 = 0.001, kb = 10000,
                     pval_threshold = 5e-8,find_proxies = TRUE,pop = "EUR",
                     harmonise_strictness = 2){
  inst_y <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                              harmonise_strictness = harmonise_strictness,
                                              find_proxies = find_proxies,
                                              pval_threshold = pval_threshold, pop = pop)
  out_y <- TwoSampleMR::extract_outcome_data(inst_y$SNP, id_outcome)
  mvdat_y <- mv_harmonise_data(inst_y, out_y)
  dat_all <- data.frame(mvdat_y$exposure_beta)
  names(dat_all) <- colnames(mvdat_y$exposure_beta)
  dat_all$y <- mvdat_y$outcome_beta
  w <- 1/mvdat_y$outcome_se^2
  if(method == "forward"){
    fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
    f1 <- stats::lm(fo1, data = dat_all, weights = w)
    fall <- stats::lm(y ~ ., data = dat_all, weights = w)
    fstep <- stats::step(f1, scope = formula(fall), direction = "forward",trace = FALSE)
    ii <- summary(fstep)$coefficients %>% rownames()
    ii <- ii[-1]
    id.select <- stringr::str_replace_all(ii, stringr::fixed("`"), "")
  }else if(method == "backward"){
    fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
    f1 <- stats::lm(fo1, data = dat_all, weights = w)
    fall <- stats::lm(y ~ ., data = dat_all, weights = w)
    fback <- stats::step(fall,scope=list(upper=formula(fall),lower=formula(f1)),
                         direction="backward",trace=FALSE)
    ii <- summary(fback)$coefficients %>% rownames()
    ii <- ii[-1]
    id.select <- stringr::str_replace_all(ii, stringr::fixed("`"), "")
  }else if(method == "both"){
    fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
    f1 <- stats::lm(fo1, data = dat_all, weights = w)
    fall <- stats::lm(y ~ ., data = dat_all, weights = w)
    fboth <- stats::step(f1,scope=list(upper=formula(fall),lower=formula(f1)),
                         direction="both",trace=FALSE)
    ii <- summary(fboth)$coefficients %>% rownames()
    ii <- ii[-1]
    id.select <- stringr::str_replace_all(ii, stringr::fixed("`"), "")
  }
  df_info[df_info$id %in% id.select,"status"] <- paste0("Select by stepwise ",method)
  return(list(id.list=id.select,trait.info=df_info))
}

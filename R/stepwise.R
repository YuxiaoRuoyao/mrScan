#' @title Use stepwise selection to do confounder selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param mvdat_y Harmonized data for the outcome produced by select_instruments step
#' @param method either forward, backward or both. Default = "forward"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import stringr
#' @export
stepwise <- function(id_exposure,id.list,df_info,mvdat_y,method = "forward"){
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
    id.select <- str_replace_all(ii, stringr::fixed("`"), "")
  }else if(method == "backward"){
    fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
    f1 <- stats::lm(fo1, data = dat_all, weights = w)
    fall <- stats::lm(y ~ ., data = dat_all, weights = w)
    fback <- stats::step(fall,scope=list(upper=formula(fall),lower=formula(f1)),
                         direction="backward",trace=FALSE)
    ii <- summary(fback)$coefficients %>% rownames()
    ii <- ii[-1]
    id.select <- str_replace_all(ii, stringr::fixed("`"), "")
  }else if(method == "both"){
    fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
    f1 <- stats::lm(fo1, data = dat_all, weights = w)
    fall <- stats::lm(y ~ ., data = dat_all, weights = w)
    fboth <- stats::step(f1,scope=list(upper=formula(fall),lower=formula(f1)),
                         direction="both",trace=FALSE)
    ii <- summary(fboth)$coefficients %>% rownames()
    ii <- ii[-1]
    id.select <- str_replace_all(ii, stringr::fixed("`"), "")
  }
  id.select <- id.select[!id.select %in% id_exposure]
  df_info[df_info$id %in% id.select,"status"] <- paste0("Select by stepwise ",method)
  return(list(id.list=id.select,trait.info=df_info))
}

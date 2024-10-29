#' @title Use stepwise selection to do confounder selection
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for zvalues, `trait_ID.p` for pvalues,
#' `trait_ID.ss` for sample sizes
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param type Input data type. It could be either "local" for local GWAS summary data
#' after LD pruning or "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param id_exposure GWAS ID of the main exposure
#' @param df_info Dataframe of trait info from previous steps
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @param method either forward, backward or both. Default = "forward"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import stringr
#' @export
stepwise <- function(dat, type, id_exposure,df_info,pval_threshold = 5e-8,
                     ss.exposure = NULL, effect_size_cutoff=0.1,
                     type_outcome = "continuous", prevalence_outcome = NULL,
                     type_exposure = NULL, prevalence_exposure = NULL,
                     method = "forward"){
  if(type == "local"){
    snp <- dat$snp
    info <- dat %>% select(snp,REF,ALT)
    beta_hat <- dat %>% select(ends_with(".beta"))
    se <- dat %>% select(ends_with(".se"))
    z <- dat %>% select(ends_with(".z"))
    p <- dat %>% select(ends_with(".p"))
    ss <- dat %>% select(ends_with(".ss"))
    af <- dat %>% select(ends_with(".af"))
    nms <- stringr::str_replace(names(z), ".z", "")
    names(beta_hat)<-names(se)<-names(z)<-names(p)<-names(ss)<-names(af)<-nms
    z.norm <- z/sqrt(ss)
    pmin <- apply(p[,-1, drop = F], 1, min)
    ix <- which(pmin < pval_threshold)
    filtered_idx <- which(rowSums(abs(data.frame(z.norm[,-1])) < effect_size_cutoff) == ncol(z.norm)-1)
    outlier_snp <- snp[-filtered_idx]
    new_ix <- intersect(ix,filtered_idx)
    filtered_SNP <- general_steiger_filtering(SNP = snp[new_ix],id.exposure = nms[-1],id.outcome = nms[1],
                                              exposure_beta = beta_hat[new_ix,-1],exposure_pval = p[new_ix,-1],
                                              exposure_se = se[new_ix,-1],outcome_beta = beta_hat[new_ix,1],
                                              outcome_pval = p[new_ix,1],outcome_se = se[new_ix,1],
                                              exposure_af = af[new_ix,-1],outcome_af = af[new_ix,1],
                                              type_outcome = type_outcome, prevalence_outcome = prevalence_outcome,
                                              type_exposure = type_exposure, prevalence_exposure = prevalence_exposure,
                                              snp_info = info[new_ix,],proxies = 0)
    final_ix <- which(snp %in% filtered_SNP)
    dat_all <- data.frame(beta_hat[final_ix,],check.names = F)
    colnames(dat_all)[1] <- "y"
    w <- 1/as.matrix(se)[final_ix,1]^2
  }
  if(type == "IEU"){
    id.exposure <- colnames(dat$exposure_beta)
    id.outcome <- dat$outname$id.outcome
    if(is.null(ss.exposure)){
      ss.exposure <- gwasinfo(id.exposure)$sample_size
    }
    z.exposure<- dat$exposure_beta/dat$exposure_se
    z.norm.exposure <- sweep(z.exposure,2,sqrt(ss.exposure),`/`)
    names(ss.exposure) <- id.exposure
    filtered_idx <- which(rowSums(abs(z.norm.exposure) < effect_size_cutoff) == ncol(z.norm.exposure))
    snp <- rownames(dat$exposure_beta)
    outlier_snp <- snp[-filtered_idx]
    filtered_SNP <- general_steiger_filtering(SNP = snp[filtered_idx],
                                              id.exposure = id.exposure,
                                              id.outcome = id.outcome,
                                              exposure_beta = dat$exposure_beta[filtered_idx,],
                                              exposure_pval = dat$exposure_pval[filtered_idx,],
                                              exposure_se = dat$exposure_se[filtered_idx,],
                                              outcome_beta = dat$outcome_beta[filtered_idx],
                                              outcome_pval = dat$outcome_pval[filtered_idx],
                                              outcome_se = dat$outcome_se[filtered_idx],
                                              type_outcome = type_outcome,
                                              prevalence_outcome = prevalence_outcome,
                                              type_exposure = type_exposure,
                                              prevalence_exposure = prevalence_exposure,
                                              proxies = 1)
    final_ix <- which(snp %in% filtered_SNP)
    dat_all <- data.frame(dat$exposure_beta[final_ix,])
    names(dat_all) <- colnames(dat$exposure_beta)
    dat_all$y <- dat$outcome_beta[final_ix]
    w <- 1/dat$outcome_se[final_ix]^2
  }
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

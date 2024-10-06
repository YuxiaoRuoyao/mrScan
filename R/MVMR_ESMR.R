#' @title Use ESMR to do MVMR analysis by locally data
#' @param dat A data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for z-values, `trait_ID.ss` for sample sizes.
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @returns A dataframe of result summary
#'
#' @import dplyr
#' @import esmr
#' @importFrom purrr map_dfc
#' @export
MVMR_ESMR <- function(dat,R_matrix,pval_threshold = 5e-8,
                      effect_size_cutoff=0.1,
                      type_outcome = "continuous", prevalence_outcome = NULL,
                      type_exposure = NULL, prevalence_exposure = NULL){
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
  o <- match(colnames(R_matrix), nms)
  beta_hat <- beta_hat[,o]
  se <- se[,o]
  z.norm <- z.norm[,o]
  i <- ncol(beta_hat)
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  filtered_idx <- which(rowSums(abs(data.frame(z.norm[,-1])) < effect_size_cutoff) == ncol(z.norm)-1)
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
  res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                            b = NA, se = NA, pvalue = NA,
                            method = paste0("ESMR_",pval_threshold))
  tryCatch({
    fit <- esmr(beta_hat_Y = beta_hat[,1],
                se_Y = se[,1],
                beta_hat_X = beta_hat[,2:i],
                se_X = se[, 2:i],
                R = R_matrix,
                variant_ix= final_ix)
    res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                              b = fit$beta$beta_m, se = fit$beta$beta_s) %>%
      mutate(pvalue = 2*pnorm(-abs(b/se)),
             method = paste0("ESMR_",pval_threshold))
    if(i <= 18){
      esmr_res <- esmr:::optimize_lpy2(fit)
      res.summary.opt <- data.frame(exposure = colnames(beta_hat)[-1],
                                    b = esmr_res$beta$beta_m, se = esmr_res$beta$beta_s) %>%
        mutate(pvalue = 2*pnorm(-abs(b/se)),
               method = paste0("ESMR_optimize_",pval_threshold))
      res.summary <- bind_rows(res.summary, res.summary.opt)
    }
  }, error = function(e){
    message("Error in MVMR_ESMR: ", e$message)
  })
  return(res.summary)
}

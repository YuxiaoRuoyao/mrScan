#' @title Use MRBEE to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for zvalues, `trait_ID.p` for pvalues,
#' `trait_ID.ss` for sample sizes
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param pleio_threshold pvalue threshold in pleiotropy detection. Default = 0
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @returns A dataframe of result summary
#'
#' @import MRBEE
#' @import dplyr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_MRBEE <- function(dat,R_matrix,pval_threshold = 5e-8,pleio_threshold = 0,type,
                       ss.exposure = NULL, effect_size_cutoff=0.1,
                       type_outcome = "continuous", prevalence_outcome = NULL,
                       type_exposure = NULL, prevalence_exposure = NULL){
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
    o <- match(colnames(R_matrix), nms)
    beta_hat <- beta_hat[,o]
    se <- se[,o]
    z.norm <- z.norm[,o]
    i <- ncol(beta_hat)
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
    # Make the last one be outcome for R matrix
    R_matrix <- R_matrix[c(nms[-1],nms[1]),c(nms[-1],nms[1])]
    if(i>2){
      fit <- MRBEE.IMRP(by=as.matrix(beta_hat)[final_ix,1],bX=as.matrix(beta_hat)[final_ix,-1],
                        byse=as.matrix(se)[final_ix,1],
                        bXse=as.matrix(se)[final_ix,-1],
                        Rxy=R_matrix,
                        pv.thres = pleio_threshold)
      res.summary <- data.frame(exposure = names(fit$theta),
                                b = fit$theta,
                                se = sqrt(diag(fit$covtheta)))
    }else{
      fit <- MRBEE.IMRP.UV(by = as.matrix(beta_hat)[final_ix,1],bx = as.matrix(beta_hat)[final_ix,-1],
                           byse = as.matrix(se)[final_ix,1],
                           bxse = as.matrix(se)[final_ix,-1],
                           Rxy=R_matrix,
                           pv.thres = pleio_threshold)
      res.summary <- data.frame(exposure = nms[-1],
                                b = fit$theta,
                                se = sqrt(fit$vartheta))
    }
    res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
    res.summary$method <- paste0("MRBEE_",pval_threshold,"_pleio_",pleio_threshold)
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
    fit <- MRBEE.IMRP(by = dat$outcome_beta[final_ix],
                      bX = dat$exposure_beta[final_ix, ],
                      byse = dat$outcome_se[final_ix],
                      bXse = dat$exposure_se[final_ix, ],
                      Rxy = diag(nrow = length(id.exposure)+1),
                      pv.thres = pleio_threshold)
    res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                              b = fit$theta,se = sqrt(diag(fit$covtheta))) %>%
      mutate(pvalue = 2*pnorm(-abs(b/se)), method = "MVMR_MRBEE")
  }
  rownames(res.summary) <- NULL
  return(list(res.summary = res.summary, outlier_SNP = outlier_snp))
}

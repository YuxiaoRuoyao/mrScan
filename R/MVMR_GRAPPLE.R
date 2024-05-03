#' @title Use GRAPPLE to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for zvalues, `trait_ID.p` for pvalues,
#' `trait_ID.ss` for sample sizes.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 1e-5
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param ss.outcome A numeric number of sample size for the outcome. You can provide it when type = "IEU". Default = NULL
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @returns A dataframe of result summary
#'
#' @import GRAPPLE
#' @import dplyr
#' @import stringr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_GRAPPLE <- function(dat,R_matrix,pval_threshold = 1e-5,type,
                         ss.exposure = NULL, ss.outcome = NULL, effect_size_cutoff=0.1,
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
    nms <- stringr::str_replace(names(z), ".z", "")
    names(beta_hat)<-names(se)<-names(z)<-names(p)<-names(ss)<-nms
    N <- apply(ss, 2, median, na.rm = TRUE)
    z.norm <- sweep(z,2,sqrt(N),`/`) %>% data.frame(check.names = F)
    se.norm <- purrr::map_dfc(ss, ~ rep(1/sqrt(.x),length.out = nrow(z))) %>%
      data.frame(check.names = F)
    o <- match(colnames(R_matrix), nms)
    z.norm <- z.norm[,o]
    se.norm <- se.norm[,o]
    i <- ncol(z.norm)
    pmin <- apply(p[,-1, drop = F], 1, min)
    ix <- which(pmin < pval_threshold)
    filtered_idx <- which(rowSums(abs(data.frame(z.norm[,-1])) < effect_size_cutoff) == ncol(z.norm)-1)
    new_ix <- intersect(ix,filtered_idx)
    filtered_SNP <- general_steiger_filtering(SNP = snp[new_ix],id.exposure = nms[-1],id.outcome = nms[1],
                                              exposure_beta = beta_hat[new_ix,-1],exposure_pval = p[new_ix,-1],
                                              exposure_se = se[new_ix,-1],outcome_beta = beta_hat[new_ix,1],
                                              outcome_pval = p[new_ix,1],outcome_se = se[new_ix,1],
                                              type_outcome = type_outcome, prevalence_outcome = prevalence_outcome,
                                              type_exposure = type_exposure, prevalence_exposure = prevalence_exposure,
                                              snp_info = info[new_ix,],proxies = 0)
    final_ix <- which(snp %in% filtered_SNP)
    grapple_dat <- data.frame(cbind(z.norm[final_ix,], se.norm[final_ix,]))
    names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(i-1)),
                            "se_out", paste0("se_exp", 1:(i-1)))
    res_and_warning <- WarningAndGrappleEst(data = grapple_dat, cor.mat = R_matrix)
    if(i > 2){
      res.summary <- data.frame(exposure=colnames(z.norm)[-1],
                                b=res_and_warning$est$beta.hat,
                                se=sqrt(diag(res_and_warning$est$beta.var)),
                                pvalue=res_and_warning$est$beta.p.value,
                                method = paste0("GRAPPLE_",pval_threshold),
                                converge = !res_and_warning$notConverge)
    }else{
      res.summary <- data.frame(exposure=colnames(z.norm)[-1],
                                b=res_and_warning$est$beta.hat,
                                se=sqrt(res_and_warning$est$beta.var),
                                pvalue=res_and_warning$est$beta.p.value,
                                method = paste0("GRAPPLE_",pval_threshold),
                                converge = !res_and_warning$notConverge)
    }
  }
  if(type == "IEU"){
    id.exposure <- colnames(dat$exposure_beta)
    id.outcome <- dat$outname$id.outcome
    if(is.null(ss.exposure)){
      ss.exposure <- gwasinfo(id.exposure)$sample_size
    }
    if(is.null(ss.outcome)){
      ss.outcome <- gwasinfo(id.outcome)$sample_size
    }
    z.exposure<- dat$exposure_beta/dat$exposure_se
    z.norm.exposure <- sweep(z.exposure,2,sqrt(ss.exposure),`/`)
    names(ss.exposure) <- id.exposure
    se.norm.exposure <- map_dfc(ss.exposure, ~ rep(1/sqrt(.x),length.out = nrow(dat$exposure_se))) %>%
      as.matrix()
    z.norm.outcome <- (dat$outcome_beta/dat$outcome_se)/sqrt(ss.outcome)
    se.norm.outcome <- rep(1/sqrt(ss.outcome),length(dat$outcome_se))
    # if(type_outcome == "binary"){
    #   ncases <- gwasinfo(id.outcome)$ncase
    #   convert_ratio <- convert_liability(k = prevalence_outcome, p = ncases/ss.outcome)
    #   z.norm.outcome <- z.norm.outcome*convert_ratio
    #   se.norm.outcome <- se.norm.outcome*convert_ratio
    # }
    # if("binary" %in% type_exposure){
    #   ix_binary <- which(type_exposure %in% "binary")
    #   ncases <- gwasinfo(id.exposure[ix_binary])$ncase
    #   convert_ratio <- convert_liability(k = prevalence_exposure[ix_binary],
    #                                      p = ncases/ss.exposure[ix_binary])
    #   if(length(ix_binary) == 1) {
    #     z.norm.exposure[, ix_binary] <- z.norm.exposure[, ix_binary]*convert_ratio
    #     se.norm.exposure[, ix_binary] <- se.norm.exposure[, ix_binary]*convert_ratio
    #   } else {
    #     z.norm.exposure[, ix_binary] <- sweep(z.norm.exposure[, ix_binary], 2, convert_ratio, `*`)
    #     se.norm.exposure[,ix_binary] <- sweep(se.norm.exposure[,ix_binary], 2, convert_ratio, `*`)
    #   }
    # }
    filtered_idx <- which(rowSums(abs(z.norm.exposure) < effect_size_cutoff) == ncol(z.norm.exposure))
    snp <- rownames(dat$exposure_beta)
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
    grapple_dat<-cbind(z.norm.exposure[final_ix,],se.norm.exposure[final_ix,],
                       z.norm.outcome[final_ix],se.norm.outcome[final_ix])
    i <- length(id.exposure)
    colnames(grapple_dat)<-c(paste0("gamma_exp",seq(1,i)),
                              paste0("se_exp",seq(1,i)),"gamma_out1","se_out1")
    res_and_warning <- WarningAndGrappleEst(data = grapple_dat)
    res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                              b = res_and_warning$est$beta.hat,
                              se = sqrt(diag(res_and_warning$est$beta.var)),
                              pvalue = res_and_warning$est$beta.p.value,
                              method = "MVMR_GRAPPLE",
                              converge = !res_and_warning$notConverge)
  }
  res.summary[which(res.summary$se > 1),"pvalue"] <- 1
  return(res.summary)
}

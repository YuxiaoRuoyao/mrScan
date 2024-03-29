#' @title Use IVW to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for zvalues, `trait_ID.p` for pvalues,
#' `trait_ID.ss` for sample sizes
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
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
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_IVW <- function(dat,pval_threshold=5e-8,type,
                     ss.exposure = NULL, ss.outcome = NULL,effect_size_cutoff = 0.1,
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
    z.norm <- sweep(z,2,sqrt(N),`/`)
    se.norm <- purrr::map_dfc(ss, ~ rep(1/sqrt(.x),length.out = nrow(z)))
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
                                              snp_info = info[new_ix,])
    final_ix <- which(snp %in% filtered_SNP)
    i <- ncol(z.norm)
    if(i>2){
      hdat <-  list(exposure_beta = as.matrix(z.norm[final_ix, 2:i,drop=FALSE]),
                    exposure_pval = as.matrix(p[final_ix, 2:i,drop = FALSE]),
                    exposure_se = as.matrix(se.norm[final_ix,2:i],drop=FALSE),
                    outcome_beta = data.frame(z.norm)[final_ix,1],
                    outcome_pval = data.frame(p)[final_ix,1],
                    outcome_se = data.frame(se.norm)[final_ix,1],
                    expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                    outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
      res_F <- mv_multiple(hdat)$result
      res_T <- mv_multiple(hdat,instrument_specific = TRUE)$result
      res_F$method <- paste0("IVW_",pval_threshold)
      res_T$method <- paste0("IVW_T_",pval_threshold)
      res <- rbind(res_F,res_T) %>% select(exposure,b,se,pval,method) %>%
        rename("pvalue" = "pval")
    }else{
      hdat <-  list(exposure_beta = as.matrix(z.norm[final_ix, 2:i,drop=FALSE]),
                    exposure_pval = as.matrix(p[final_ix, 2:i,drop = FALSE]),
                    exposure_se = as.matrix(se.norm[final_ix,2:i],drop=FALSE),
                    outcome_beta = data.frame(z.norm)[final_ix,1],
                    outcome_pval = data.frame(p)[final_ix,1],
                    outcome_se = data.frame(se.norm)[final_ix,1],
                    expname = data.frame(id.exposure = nms[-1], exposure = nms[-1]),
                    outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
      res_F <- mv_multiple(hdat)$result
      res_F$method <- paste0("IVW_",pval_threshold)
      res <- res_F %>% select(exposure,b,se,pval,method) %>%
        rename("pvalue" = "pval")
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
    if(type_outcome == "binary"){
      ncases <- gwasinfo(id.outcome)$ncase
      convert_ratio <- convert_liability(k = prevalence_outcome, p = ncases/ss.outcome)
      z.norm.outcome <- z.norm.outcome*convert_ratio
      se.norm.outcome <- se.norm.outcome*convert_ratio
    }
    if("binary" %in% type_exposure){
      ix_binary <- which(type_exposure %in% "binary")
      ncases <- gwasinfo(id.exposure[ix_binary])$ncase
      convert_ratio <- convert_liability(k = prevalence_exposure[ix_binary],
                                         p = ncases/ss.exposure[ix_binary])
      if(length(ix_binary) == 1) {
        z.norm.exposure[, ix_binary] <- z.norm.exposure[, ix_binary]*convert_ratio
        se.norm.exposure[, ix_binary] <- se.norm.exposure[, ix_binary]*convert_ratio
      } else {
        z.norm.exposure[, ix_binary] <- sweep(z.norm.exposure[, ix_binary], 2, convert_ratio, `*`)
        se.norm.exposure[,ix_binary] <- sweep(se.norm.exposure[,ix_binary], 2, convert_ratio, `*`)
      }
    }
    filtered_idx <- which(rowSums(abs(z.norm.exposure) < effect_size_cutoff) == ncol(z.norm.exposure))

    dat_norm <- list(exposure_beta = z.norm.exposure[filtered_idx, ],
                     exposure_pval = dat$exposure_pval[filtered_idx, ],
                     exposure_se = se.norm.exposure[filtered_idx, ],
                     outcome_beta = z.norm.outcome[filtered_idx],
                     outcome_pval = dat$outcome_pval[filtered_idx],
                     outcome_se = se.norm.outcome[filtered_idx],
                     expname = dat$expname, outname = dat$outname)
    res <- mv_multiple(dat_norm)$result %>% select(id.exposure,id.outcome,b,se,pval) %>%
      rename("pvalue" = "pval") %>% mutate(method = "MVMR_IVW")
  }
  return(res)
}

#' @title Calculates the conditional F-statistic to assess instrument strength and filter traits
#' @param dat A data frame of combined GWAS summary data after LD pruning.
#' For local data, the required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' For IEU data, it should be output from TwoSampleMR::mv_harmonise_data().
#' @param dat_type Either "local" or "IEU". Default  = "local"
#' @param R_matrix Pairwise sample overlap matrix among traits. Default = NULL
#' @param df_info Dataframe of trait info from previous step
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param F_threshold F-statistic cutoff. Default = 5
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param min_instruments minimum number of instruments. Default = 3
#' @param Filter Whether perform trait filtering based on F-stats and F_threshold. Default = FALSE
#' @param extra_traits trait ID you want to include no matter the instrument strength. Default = "None"
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @param ss.exposure A vector of sample size for exposures. You can provide it when dat_type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param df_af_out A dataframe of allele frequency of the outcome.
#' It contains columns `SNP`,`eaf.outcome`,`beta.outcome`,`id.outcome`. Default = NULL
#' @param df_af_exp A list for allele frequency matrix for each exposure.
#' Each dataframe contains columns `SNP`,`eaf.exposure`,`beta.exposure`,`id.exposure`.
#' Each element in the list is the dataframe for each exposure. Default = NULL
#' @returns A list of selected traits, a dataframe of conditional instrument strength and a dataframe of trait info
#'
#' @import dplyr
#' @import MVMR
#' @import ieugwasr
#' @importFrom purrr map_dfr
#' @export
strength_filter <- function(dat,dat_type = "local",R_matrix = NULL,df_info,
                            pval_threshold = 5e-8,F_threshold = 5,
                            effect_size_cutoff = 0.1,
                            min_instruments = 3,
                            Filter = FALSE, extra_traits = "None",
                            type_outcome = "continuous", prevalence_outcome = NULL,
                            type_exposure = NULL, prevalence_exposure = NULL,
                            ss.exposure = NULL, df_af_out = NULL, df_af_exp = NULL){
  if(dat_type == "local"){
    snp <- dat$snp
    info <- dat %>% select(snp,REF,ALT)
    beta_hat <- dat %>% select(ends_with(".beta"))
    se <- dat %>% select(ends_with(".se"))
    z <- dat %>% select(ends_with(".z"))
    p <- dat %>% select(ends_with(".p"))
    ss <- dat %>% select(ends_with(".ss"))
    af <- dat %>% select(ends_with(".af"))
    nms <- stringr::str_replace(names(beta_hat), ".beta", "")
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
    exp_se <- as.matrix(se[final_ix,2:i])
    if(i > 2){
      F.data<- format_mvmr(BXGs = as.matrix(beta_hat[final_ix, 2:i]),
                           BYG = beta_hat[final_ix,1],
                           seBXGs = exp_se,
                           seBYG = se[final_ix,1],
                           RSID = snp[final_ix])
      sigmalist <- vector("list", length(final_ix))
      if(!is.null(R_matrix)){
        omega <- R_matrix[2:i,2:i]
        for (j in 1:length(final_ix)) {
          se_matrix <- diag(exp_se[j,],nrow = i-1)
          sigmalist[[j]] <- se_matrix %*% omega %*% se_matrix
        }
        sres <- data.frame(t(strength_mvmr(r_input = F.data, gencov = sigmalist)))
      }else{
        sres <- data.frame(t(strength_mvmr(r_input = F.data, gencov = 0)))
      }
    }else{
      MRInputObject <- MendelianRandomization::mr_input(bx = as.matrix(beta_hat)[final_ix, 2:i],
                                                        by = as.matrix(beta_hat)[final_ix,1],
                                                        bxse = as.matrix(se)[final_ix,2:i],
                                                        byse = as.matrix(se)[final_ix,1])
      IVWObject <- MendelianRandomization::mr_ivw(MRInputObject)
      sres <- data.frame(F.statistic = IVWObject@Fstat)
    }
    sres$id<-colnames(beta_hat[-1])
  }
  if(dat_type == "IEU"){
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
    info_outcome <- ieugwasr::gwasinfo(id.outcome)
    if (nrow(info_outcome) == 0) {
      outcome_af <- data.frame(df_af_out[filtered_idx, "eaf.outcome"])
      ncase_outcome <- unique(df_af_out$ncase.outcome)
      ncontrol_outcome <- unique(df_af_out$ncontrol.outcome)
      samplesize_outcome <- unique(df_af_out$samplesize.outcome)
    } else {
      outcome_af <- NULL
      ncase_outcome <- NULL
      ncontrol_outcome <- NULL
      samplesize_outcome <- NULL
    }
    exposure_af <- data.frame(do.call(cbind, lapply(df_af_exp, function(df) {
      df$eaf.exposure[filtered_idx]
    })))
    colnames(exposure_af) <- names(df_af_exp)
    ncase_exposure <- c()
    ncontrol_exposure <- c()
    samplesize_exposure <- c()
    for (i in seq_along(id.exposure)) {
      info_exposure <- ieugwasr::gwasinfo(id.exposure[i])
      if (nrow(info_exposure) == 0) {
        ncase_exposure[i] <- unique(df_af_exp[[i]]$ncase.exposure)
        ncontrol_exposure[i] <- unique(df_af_exp[[i]]$ncontrol.exposure)
        samplesize_exposure[i] <- unique(df_af_exp[[i]]$samplesize.exposure)
      } else {
        ncase_exposure[i] <- NA
        ncontrol_exposure[i] <- NA
        samplesize_exposure[i] <- NA
      }
    }
    filtered_SNP <- general_steiger_filtering(SNP = snp[filtered_idx],
                                              id.exposure = id.exposure,
                                              id.outcome = id.outcome,
                                              exposure_beta = dat$exposure_beta[filtered_idx,],
                                              exposure_pval = dat$exposure_pval[filtered_idx,],
                                              exposure_se = dat$exposure_se[filtered_idx,],
                                              exposure_af = exposure_af,
                                              outcome_beta = dat$outcome_beta[filtered_idx],
                                              outcome_pval = dat$outcome_pval[filtered_idx],
                                              outcome_se = dat$outcome_se[filtered_idx],
                                              outcome_af = outcome_af,
                                              type_outcome = type_outcome,
                                              prevalence_outcome = prevalence_outcome,
                                              type_exposure = type_exposure,
                                              prevalence_exposure = prevalence_exposure,
                                              proxies = 1,
                                              ncase_outcome = ncase_outcome,
                                              ncontrol_outcome = ncontrol_outcome,
                                              samplesize_outcome = samplesize_outcome,
                                              ncase_exposure = ncase_exposure,
                                              ncontrol_exposure = ncontrol_exposure,
                                              samplesize_exposure = samplesize_exposure)
    final_ix <- which(snp %in% filtered_SNP)
    if(length(final_ix) > min_instruments){
      F.data <- format_mvmr(BXGs = as.matrix(dat$exposure_beta[final_ix,]),
                            BYG = dat$outcome_beta[final_ix],
                            seBXGs = as.matrix(dat$exposure_se[final_ix,]),
                            seBYG = dat$outcome_se[final_ix],
                            RSID = rownames(dat$exposure_beta)[final_ix])
      sres <- data.frame(t(strength_mvmr(r_input = F.data, gencov = 0)))
    }
    else{
      sres <- data.frame(F.statistic = rep(NA, ncol(dat$exposure_beta)))
      print("Not enough instruments!")
    }
    sres$id <- colnames(dat$exposure_beta)
  }
  if(Filter == TRUE){
    select.id <- sres %>% filter(F.statistic > F_threshold) %>% pull(id)
    if(extra_traits != "None"){
      select.id <- unique(c(select.id,extra_traits))
    }
    other.id <- sres$id[!sres$id %in% select.id]
    df_info[df_info$id %in% other.id,"status"] <- "Delete due to weak instruments strength"
    df_strength <- data.frame(sres) %>% left_join(df_info[,c("id","trait")])
    return(list(id.list = select.id,df_strength = df_strength,trait.info=df_info))
  }else{
    df_strength <- data.frame(sres) %>% left_join(df_info[,c("id","trait")])
    return(df_strength)
  }
}

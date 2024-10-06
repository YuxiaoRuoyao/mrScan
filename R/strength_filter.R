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
#' @param Filter Whether perform trait filtering based on F-stats and F_threshold. Default = FALSE
#' @param extra_traits trait ID you want to include no matter the instrument strength. Default = "None"
#' @param type_outcome It could be either "continuous" or "binary". Default = "continuous"
#' @param prevalence_outcome Outcome prevalence. It should be input if the outcome is a binary trait. Default = NULL
#' @param type_exposure A vector for the type of exposures. The order should be exactly matched
#' with exposures. eg. c("continuous","binary","continuous") for the second exposure is a binary trait
#' @param prevalence_exposure A vector for prevalence of exposures. The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(NA, 0.1, NA)
#' @returns A list of selected traits, a dataframe of conditional instrument strength and a dataframe of trait info
#'
#' @import dplyr
#' @import MVMR
#' @importFrom purrr map_dfr
#' @export
strength_filter <- function(dat,dat_type = "local",R_matrix = NULL,df_info,
                            pval_threshold = 5e-8,F_threshold = 5,
                            effect_size_cutoff = 0.1,
                            Filter = FALSE, extra_traits = "None",
                            type_outcome = "continuous", prevalence_outcome = NULL,
                            type_exposure = NULL, prevalence_exposure = NULL){
  if(dat_type == "local"){
    snp <- data.frame(dat$snp)
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
    beta_hat <- data.frame(beta_hat[, o],check.names = F)
    se <- data.frame(se[, o],check.names = F)
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
    F.data<- format_mvmr(BXGs = as.matrix(beta_hat[final_ix, 2:i]),
                         BYG = beta_hat[final_ix,1],
                         seBXGs = exp_se,
                         seBYG = se[final_ix,1],
                         RSID = snp[final_ix,1])
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
    F.data <- format_mvmr(BXGs = as.matrix(dat$exposure_beta[final_ix,]),
                         BYG = dat$outcome_beta[final_ix],
                         seBXGs = as.matrix(dat$exposure_se[final_ix,]),
                         seBYG = dat$outcome_se[final_ix],
                         RSID = rownames(dat$exposure_beta)[final_ix])
    sres <- data.frame(t(strength_mvmr(r_input = F.data, gencov = 0)))
    sres$id<-colnames(dat$exposure_beta)
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

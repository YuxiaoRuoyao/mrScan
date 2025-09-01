#' @title Conduct bidirection MR between traits
#' @param ex_dat1 Dataframe of instruments for the main exposure or outcome. Output from TwoSampleMR::extract_instruments()
#' @param ex_dat2 Dataframe of instruments for a candidate confounder trait. Output from TwoSampleMR::extract_instruments()
#' @param min_instruments minimum number of instruments for candidate traits. Default = 3
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param R2_cutoff R2 cutoff for duplicated traits with X or Y. Default = 0.85
#' @param type_list A vector for the type of traits. The order should be exactly matched
#' with `ID1`, `ID2`. ID1 is the exposure of ex_dat1 and ID2 is the exposure of ex_dat2.
#' eg. c("binary","continuous") means that ID1 is a binary trait.
#' @param prevalence_list A list for prevalence of traits. The order should
#' be exactly matched with `ID1`, `ID2`. eg. list(0.1, NULL) means the prevalence of
#' ID1 is 0.1. For continuous trait, just write NULL.
#' @param df Optional local file. Default = NULL
#' @returns A list contain bidirection estimates and traits correlation
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @export
bidirection_mr <- function(ex_dat1,ex_dat2,min_instruments=3,effect_size_cutoff=0.1,R2_cutoff=0.85,
                           type_list = c("continuous","continuous"), prevalence_list = NULL,
                           df = NULL){
  ID1 <- unique(ex_dat1$id.exposure)
  ID2 <- unique(ex_dat2$id.exposure)
  if(is.null(ex_dat2)){
    return(NULL)
  }else if(!is.null(ex_dat2) & nrow(ex_dat2) < min_instruments){
    return(NULL)
  }else{
    out_dat1 <- extract_outcome_data(snps = ex_dat1$SNP,outcomes = ID2)
    info_ID1 <- ieugwasr::gwasinfo(ID1)
    if(nrow(info_ID1) != 0){
      out_dat2 <- extract_outcome_data(snps = ex_dat2$SNP,outcomes = ID1)
    }else{
      out_dat2 <- format_data(as.data.frame(df), type = "outcome",
                              snps = ex_dat2$SNP, snp_col = "hm_rsid",
                              beta_col = "hm_beta", se_col = "standard_error",
                              eaf_col = "hm_effect_allele_frequency",
                              effect_allele_col = "hm_effect_allele", other_allele_col = "hm_other_allele",
                              pval_col = "p_value", ncase_col = "num_cases",
                              ncontrol_col = "num_controls", samplesize_col = "sample_size",
                              chr_col = "hm_chrom", pos_col = "hm_pos")
    }
    dat_1_2 <- harmonise_data(ex_dat1, out_dat1) %>% filter(mr_keep == TRUE)
    dat_2_1 <- harmonise_data(ex_dat2, out_dat2) %>% filter(mr_keep == TRUE)
    if(sum(is.na(dat_1_2$samplesize.exposure)) != 0){
      dat_1_2$samplesize.exposure <- gwasinfo(ID1)$sample_size
    }else if(sum(is.na(dat_2_1$samplesize.exposure)) != 0){
      dat_2_1$samplesize.exposure <- gwasinfo(ID2)$sample_size
    }
    X_1_2 <- dat_1_2 %>%
      rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      select(SNP, beta1, beta2)
    X_2_1 <- dat_2_1 %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    cor_vals<- data.frame(id1 = ID1,id2 = ID2, cor = with(X, cor(beta1, beta2)))
    if(abs(cor_vals$cor) > R2_cutoff){
      return(list(mr12 = NULL, mr21 = NULL, cor = cor_vals))
    }else{
      dat_1_2 <- dat_1_2 %>%
        mutate(z.norm.exposure = (beta.exposure/se.exposure)/sqrt(samplesize.exposure),
               se.norm.exposure = 1/sqrt(samplesize.exposure))
      filtered_idx_1_2 <- which(abs(dat_1_2$z.norm.exposure) < effect_size_cutoff)
      dat_2_1 <- dat_2_1 %>%
        mutate(z.norm.exposure = (beta.exposure/se.exposure)/sqrt(samplesize.exposure),
               se.norm.exposure = 1/sqrt(samplesize.exposure))
      filtered_idx_2_1 <- which(abs(dat_2_1$z.norm.exposure) < effect_size_cutoff)
      outlier_snp_1_2 <- dat_1_2$SNP[-filtered_idx_1_2]
      outlier_snp_2_1 <- dat_2_1$SNP[-filtered_idx_2_1]
      filtered_SNP_1_2 <- general_steiger_filtering(SNP = dat_1_2$SNP[filtered_idx_1_2],
                                                    id.exposure = ID1,id.outcome = ID2,
                                                    exposure_beta = data.frame(dat_1_2$beta.exposure[filtered_idx_1_2]),
                                                    exposure_pval = data.frame(dat_1_2$pval.exposure[filtered_idx_1_2]),
                                                    exposure_se = data.frame(dat_1_2$se.exposure[filtered_idx_1_2]),
                                                    exposure_af = data.frame(dat_1_2$eaf.exposure[filtered_idx_1_2]),
                                                    outcome_beta = dat_1_2$beta.outcome[filtered_idx_1_2],
                                                    outcome_pval = dat_1_2$pval.outcome[filtered_idx_1_2],
                                                    outcome_se = dat_1_2$se.outcome[filtered_idx_1_2],
                                                    outcome_af = data.frame(dat_1_2$eaf.outcome[filtered_idx_1_2]),
                                                    type_outcome = type_list[2],
                                                    prevalence_outcome = prevalence_list[[2]],
                                                    type_exposure = type_list[1],
                                                    prevalence_exposure = prevalence_list[[1]],
                                                    ncase_exposure = unique(dat_1_2$ncase.exposure),
                                                    ncontrol_exposure = unique(dat_1_2$ncontrol.exposure),
                                                    samplesize_exposure = unique(dat_1_2$samplesize.exposure),
                                                    ncase_outcome = unique(dat_1_2$ncase.outcome),
                                                    ncontrol_outcome = unique(dat_1_2$ncontrol.outcome),
                                                    samplesize_outcome = unique(dat_1_2$samplesize.outcome),
                                                    proxies = 1)
      filtered_SNP_2_1 <- general_steiger_filtering(SNP = dat_2_1$SNP[filtered_idx_2_1],
                                                    id.exposure = ID2,id.outcome = ID1,
                                                    exposure_beta = data.frame(dat_2_1$beta.exposure[filtered_idx_2_1]),
                                                    exposure_pval = data.frame(dat_2_1$pval.exposure[filtered_idx_2_1]),
                                                    exposure_se = data.frame(dat_2_1$se.exposure[filtered_idx_2_1]),
                                                    exposure_af = data.frame(dat_2_1$eaf.exposure[filtered_idx_2_1]),
                                                    outcome_beta = dat_2_1$beta.outcome[filtered_idx_2_1],
                                                    outcome_pval = dat_2_1$pval.outcome[filtered_idx_2_1],
                                                    outcome_se = dat_2_1$se.outcome[filtered_idx_2_1],
                                                    outcome_af = data.frame(dat_2_1$eaf.outcome[filtered_idx_2_1]),
                                                    type_outcome = type_list[1],
                                                    prevalence_outcome = prevalence_list[[1]],
                                                    type_exposure = type_list[2],
                                                    prevalence_exposure = prevalence_list[[2]],
                                                    ncase_exposure = unique(dat_2_1$ncase.exposure),
                                                    ncontrol_exposure = unique(dat_2_1$ncontrol.exposure),
                                                    samplesize_exposure = unique(dat_2_1$samplesize.exposure),
                                                    ncase_outcome = unique(dat_2_1$ncase.outcome),
                                                    ncontrol_outcome = unique(dat_2_1$ncontrol.outcome),
                                                    samplesize_outcome = unique(dat_2_1$samplesize.outcome),
                                                    proxies = 1)
      final_ix_1_2 <- which(dat_1_2$SNP %in% filtered_SNP_1_2)
      final_ix_2_1 <- which(dat_2_1$SNP %in% filtered_SNP_2_1)
      methods <- list(MR_IVW = MR_IVW, MR_GRAPPLE = MR_GRAPPLE, MR_MRBEE = MR_MRBEE)
      if(length(final_ix_1_2) < min_instruments){
        res_1_2 <- NULL
      }else{
          params1 <- list(id.exposure = ID1,id.outcome = ID2,
                          beta.exposure = dat_1_2$beta.exposure[final_ix_1_2],
                          beta.outcome = dat_1_2$beta.outcome[final_ix_1_2],
                          se.exposure = dat_1_2$se.exposure[final_ix_1_2],
                          se.outcome = dat_1_2$se.outcome[final_ix_1_2])
          res_1_2 <- lapply(methods, function(f, params) {
            do.call(f, params)}, params = params1) %>% bind_rows()
      }
      if(length(final_ix_2_1) < min_instruments){
        res_2_1 <- NULL
      } else {
        params2 <- list(id.exposure = ID2,id.outcome = ID1,
                        beta.exposure = dat_2_1$beta.exposure[final_ix_2_1],
                        beta.outcome = dat_2_1$beta.outcome[final_ix_2_1],
                        se.exposure = dat_2_1$se.exposure[final_ix_2_1],
                        se.outcome = dat_2_1$se.outcome[final_ix_2_1])
        res_2_1 <- lapply(methods, function(f, params) {
          do.call(f, params)}, params = params2) %>% bind_rows()
      }
      return(list(mr12 = res_1_2, mr21 = res_2_1, cor = cor_vals,
                  mr12_outlier_SNP = outlier_snp_1_2, mr21_outlier_SNP = outlier_snp_2_1))
    }
  }
}

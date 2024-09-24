#' @title Conduct bidirection MR between traits
#' @param ex_dat1 Dataframe of instruments for the main exposure or outcome. Output from TwoSampleMR::extract_instruments()
#' @param ex_dat2 Dataframe of instruments for a candidate confounder trait. Output from TwoSampleMR::extract_instruments()
#' @param min_instruments minimum number of instruments for candidate traits. Default = 3
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.05
#' @param R2_cutoff R2 cutoff for duplicated traits with X or Y. Default = 0.85
#' @param type_list A vector for the type of traits. The order should be exactly matched
#' with `ID1`, `ID2`. ID1 is the exposure of ex_dat1 and ID2 is the exposure of ex_dat2.
#' eg. c("binary","continuous") means that ID1 is a binary trait.
#' @param prevalence_list A list for prevalence of traits. The order should
#' be exactly matched with `ID1`, `ID2`. For continuous trait, just write NULL. eg. list(0.1, NULL)
#' @returns A list contain bidirection estimates and traits correlation
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @export
bidirection_mr <- function(ex_dat1,ex_dat2,min_instruments=3,effect_size_cutoff=0.1,R2_cutoff=0.85,
                           type_list = c("continuous","continuous"), prevalence_list = NULL){
  ID1 <- unique(ex_dat1$id.exposure)
  ID2 <- unique(ex_dat2$id.exposure)
  if(is.null(ex_dat2)){
    return(NULL)
  }else if(!is.null(ex_dat2) & nrow(ex_dat2) < min_instruments){
    return(NULL)
  }else{
    out_dat1 <- extract_outcome_data(snps = ex_dat1$SNP,outcomes = ID2)
    out_dat2 <- extract_outcome_data(snps = ex_dat2$SNP,outcomes = ID1)
    dat_1_2 <- harmonise_data(ex_dat1, out_dat1) %>% filter(mr_keep == TRUE)
    dat_2_1 <- harmonise_data(ex_dat2, out_dat2) %>% filter(mr_keep == TRUE)
    if(sum(is.na(dat_1_2$samplesize.exposure)) != 0){
      dat_1_2$samplesize.exposure <- dat_2_1$samplesize.outcome <- gwasinfo(ID1)$sample_size
    }else if(sum(is.na(dat_1_2$samplesize.outcome)) != 0){
      dat_1_2$samplesize.outcome <- dat_2_1$samplesize.exposure <- gwasinfo(ID2)$sample_size
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
               se.norm.exposure = 1/sqrt(samplesize.exposure),
               z.norm.outcome = (beta.outcome/se.outcome)/sqrt(samplesize.outcome),
               se.norm.outcome = 1/sqrt(samplesize.outcome)) %>%
        filter(abs(z.norm.exposure) < effect_size_cutoff)
      dat_2_1 <- dat_2_1 %>%
        mutate(z.norm.exposure = (beta.exposure/se.exposure)/sqrt(samplesize.exposure),
               se.norm.exposure = 1/sqrt(samplesize.exposure),
               z.norm.outcome = (beta.outcome/se.outcome)/sqrt(samplesize.outcome),
               se.norm.outcome = 1/sqrt(samplesize.outcome)) %>%
        filter(abs(z.norm.exposure) < effect_size_cutoff)
      filtered_SNP_1_2 <- general_steiger_filtering(SNP = dat_1_2$SNP,
                                                    id.exposure = ID1,id.outcome = ID2,
                                                    exposure_beta = data.frame(dat_1_2$beta.exposure),
                                                    exposure_pval = data.frame(dat_1_2$pval.exposure),
                                                    exposure_se = data.frame(dat_1_2$se.exposure),
                                                    exposure_af = data.frame(dat_1_2$eaf.exposure),
                                                    outcome_beta = dat_1_2$beta.outcome,
                                                    outcome_pval = dat_1_2$pval.outcome,
                                                    outcome_se = dat_1_2$se.outcome,
                                                    outcome_af = data.frame(dat_1_2$eaf.outcome),
                                                    type_outcome = type_list[2],
                                                    prevalence_outcome = prevalence_list[[2]],
                                                    type_exposure = type_list[1],
                                                    prevalence_exposure = prevalence_list[[1]],
                                                    proxies = 1)
      filtered_SNP_2_1 <- general_steiger_filtering(SNP = dat_2_1$SNP,
                                                    id.exposure = ID2,id.outcome = ID1,
                                                    exposure_beta = data.frame(dat_2_1$beta.exposure),
                                                    exposure_pval = data.frame(dat_2_1$pval.exposure),
                                                    exposure_se = data.frame(dat_2_1$se.exposure),
                                                    exposure_af = data.frame(dat_2_1$eaf.exposure),
                                                    outcome_beta = dat_2_1$beta.outcome,
                                                    outcome_pval = dat_2_1$pval.outcome,
                                                    outcome_se = dat_2_1$se.outcome,
                                                    outcome_af = data.frame(dat_2_1$eaf.outcome),
                                                    type_outcome = type_list[1],
                                                    prevalence_outcome = prevalence_list[[1]],
                                                    type_exposure = type_list[2],
                                                    prevalence_exposure = prevalence_list[[2]],
                                                    proxies = 1)
      dat_1_2 <- dat_1_2 %>% filter(SNP %in% filtered_SNP_1_2)
      dat_2_1 <- dat_2_1 %>% filter(SNP %in% filtered_SNP_2_1)
      methods <- list(MR_IVW = MR_IVW, MR_GRAPPLE = MR_GRAPPLE, MR_MRBEE = MR_MRBEE)
      if(nrow(dat_1_2) < min_instruments){
        res_1_2 <- NULL
      }else{
          params1 <- list(id.exposure = ID1,id.outcome = ID2,
                          z.norm.exposure = dat_1_2$z.norm.exposure,
                          z.norm.outcome = dat_1_2$z.norm.outcome,
                          se.norm.exposure = dat_1_2$se.norm.exposure,
                          se.norm.outcome = dat_1_2$se.norm.outcome)
          res_1_2 <- lapply(methods, function(f, params) {
            do.call(f, params)}, params = params1) %>% bind_rows()
      }
      if(nrow(dat_2_1) < min_instruments){
        res_2_1 <- NULL
      } else {
        params2 <- list(id.exposure = ID2,id.outcome = ID1,
                        z.norm.exposure = dat_2_1$z.norm.exposure,
                        z.norm.outcome = dat_2_1$z.norm.outcome,
                        se.norm.exposure = dat_2_1$se.norm.exposure,
                        se.norm.outcome = dat_2_1$se.norm.outcome)
        res_2_1 <- lapply(methods, function(f, params) {
          do.call(f, params)}, params = params2) %>% bind_rows()
      }
      if (is.null(res_1_2) && is.null(res_2_1)) {
        return(NULL)
      }
      return(list(mr12 = res_1_2, mr21 = res_2_1, cor = cor_vals))
    }
  }
}

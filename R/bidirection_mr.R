#' @title Conduct bidirection MR between traits
#' @param ex_dat1 Dataframe of instruments for the main exposure or outcome. Output from TwoSampleMR::extract_instruments()
#' @param ex_dat2 Dataframe of instruments for a candidate confounder trait. Output from TwoSampleMR::extract_instruments()
#' @param min_instruments minimum number of instruments for candidate traits. Default = 3
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.05
#' @param R2_cutoff R2 cutoff for duplicated traits with X or Y. Default = 0.85
#' @returns A list contain bidirection estimates and traits correlation
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import ieugwasr
#' @export
bidirection_mr <- function(ex_dat1,ex_dat2,min_instruments=3,effect_size_cutoff=0.05,R2_cutoff=0.85){
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
      if(nrow(dat_1_2) < min_instruments | nrow(dat_2_1) < min_instruments){
        return(NULL)
      }else{
        dat_1_2 <- dat_1_2 %>% steiger_filtering() %>% filter(steiger_dir == TRUE)
        dat_2_1 <- dat_2_1 %>% steiger_filtering() %>% filter(steiger_dir == TRUE)
        if(nrow(dat_1_2) < min_instruments | nrow(dat_2_1) < min_instruments){
          return(NULL)
        }else{
          methods <- list(MR_IVW = MR_IVW, MR_GRAPPLE = MR_GRAPPLE, MR_MRBEE = MR_MRBEE)
          params1 <- list(id.exposure = unique(dat_1_2$id.exposure),
                          id.outcome = unique(dat_1_2$id.outcome),
                          z.norm.exposure = dat_1_2$z.norm.exposure,
                          z.norm.outcome = dat_1_2$z.norm.outcome,
                          se.norm.exposure = dat_1_2$se.norm.exposure,
                          se.norm.outcome = dat_1_2$se.norm.outcome)
          params2 <- list(id.exposure = unique(dat_2_1$id.exposure),
                          id.outcome = unique(dat_2_1$id.outcome),
                          z.norm.exposure = dat_2_1$z.norm.exposure,
                          z.norm.outcome = dat_2_1$z.norm.outcome,
                          se.norm.exposure = dat_2_1$se.norm.exposure,
                          se.norm.outcome = dat_2_1$se.norm.outcome)
          res_1_2 <- lapply(methods, function(f, params) {
            do.call(f, params)}, params = params1) %>% bind_rows()
          res_2_1 <- lapply(methods, function(f, params) {
            do.call(f, params)}, params = params2) %>% bind_rows()
          return(list(mr12 = res_1_2, mr21 = res_2_1, cor = cor_vals))
        }
      }
    }
  }
}

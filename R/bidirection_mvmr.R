#' @title Conduct bidirection MVMR between traits
#' @param ex_dat1 Dataframe of instruments for the main exposure or outcome and extra trait (X/Y + M). Output from TwoSampleMR::mv_extract_exposures()
#' @param ex_dat2 Dataframe of instruments for a candidate confounder trait (Z). Output from TwoSampleMR::extract_instruments()
#' @param ex_dat3 Dataframe of instruments for a candidate confounder trait and extra trait (Z + M). Output from TwoSampleMR::mv_extract_exposures()
#' @param ex_dat4 Dataframe of instruments for the main exposure or outcome (X/Y). Output from TwoSampleMR::extract_instruments()
#' @param min_instruments minimum number of instruments for candidate traits. Default = 3
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @param R2_cutoff R2 cutoff for duplicated traits with X or Y. Default = 0.85
#' @param df_info Dataframe of trait info containing sample sizes. The required columns include
#' `id` for trait ID, `sample_size` for sample sizes. Default = NULL
#' @param type_list A vector for the type of traits (X/Y + M). The order should be exactly matched
#' with traits. eg. c("binary","continuous") for the first exposure is a binary trait
#' @param prevalence_list A vector for prevalence of traits (X/Y + M). The order should
#' be exactly matched with exposures. For continuous trait, just write NA. eg. c(0.1, NA)
#' @returns A list contain bidirection estimates and traits correlation
#'
#' @import TwoSampleMR
#' @import dplyr
#' @import GRAPPLE
#' @export
bidirection_mvmr <- function(ex_dat1,ex_dat2,ex_dat3,ex_dat4,min_instruments = 3,
                             effect_size_cutoff = 0.1,R2_cutoff=0.85,df_info = NULL,
                             type_list = c("continuous","continuous"), prevalence_list = NULL){
  ID1 <- unique(ex_dat4$id.exposure) # X/Y
  ID2 <- unique(ex_dat2$id.exposure) # Z
  ID3 <- unique(ex_dat1$id.exposure)[-1] # M
  if(is.null(ex_dat2)){
    return(NULL)
  }else if(!is.null(ex_dat2) & nrow(ex_dat2) < min_instruments){
    return(NULL)
  }else{
    out_dat1 <- extract_outcome_data(snps = ex_dat4$SNP,outcomes = ID2)
    out_dat2 <- extract_outcome_data(snps = ex_dat2$SNP,outcomes = ID1)
    dat_1_2 <- harmonise_data(ex_dat4, out_dat1)
    dat_2_1 <- harmonise_data(ex_dat2, out_dat2)
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
      out_1_2 <- extract_outcome_data(snps = ex_dat1$SNP,outcomes = ID2)
      out_3_4 <- extract_outcome_data(snps = ex_dat3$SNP,outcomes = ID1)
      mvdat_1 <- mv_harmonise_data(ex_dat1,out_1_2) # ID1 + ID3 to ID2
      mvdat_2 <- mv_harmonise_data(ex_dat3,out_3_4) # ID2 + ID3 to ID1
      if(!is.null(df_info) & sum(c(ID1, ID2, ID3) %in% df_info$id) == 3){
        ss <- df_info %>% filter(id %in% c(ID1, ID2, ID3)) %>%
          arrange(match(id, c(ID1, ID2, ID3))) %>% pull(sample_size)
        params1 <- list(dat = mvdat_1, type = "IEU", ss.exposure = ss[-2],
                        effect_size_cutoff = effect_size_cutoff,
                        type_exposure = type_list, prevalence_exposure = prevalence_list)
        params2 <- list(dat = mvdat_2, type = "IEU", ss.exposure = ss[-1],
                        effect_size_cutoff = effect_size_cutoff,
                        type_outcome = type_list[1], prevalence_outcome = prevalence_list[1])
      }else{
        params1 <- list(dat = mvdat_1, type = "IEU",effect_size_cutoff = effect_size_cutoff,
                        type_exposure = type_list, prevalence_exposure = prevalence_list)
        params2 <- list(dat = mvdat_2, type = "IEU",effect_size_cutoff = effect_size_cutoff,
                        type_outcome = type_list[1], prevalence_outcome = prevalence_list[1])
      }
      methods <- list(MVMR_IVW = MVMR_IVW, MVMR_GRAPPLE = MVMR_GRAPPLE,
                      MVMR_MRBEE = MVMR_MRBEE)
      mr12_results <- lapply(methods, function(f, params) {
        res <- do.call(f, params)
        return(list(summary = res[[1]], outlier_SNP = res[[2]]))
      }, params = params1)
      mr21_results <- lapply(methods, function(f, params) {
        res <- do.call(f, params)
        return(list(summary = res[[1]], outlier_SNP = res[[2]]))
      }, params = params2)
      res1 <- lapply(mr12_results, `[[`, "summary") %>% bind_rows() %>% filter(id.exposure == ID1)
      mr12_outliers <- lapply(mr12_results, `[[`, "outlier_SNP") %>% unlist() %>% unique()
      res2 <- lapply(mr21_results, `[[`, "summary") %>% bind_rows() %>% filter(id.exposure == ID2)
      mr21_outliers <- lapply(mr21_results, `[[`, "outlier_SNP") %>% unlist() %>% unique()
      return(list(mr12 = res1, mr21 = res2, cor = cor_vals,
                  mr12_outlier_SNP = mr12_outliers, mr21_outlier_SNP = mr21_outliers))
    }
  }
}

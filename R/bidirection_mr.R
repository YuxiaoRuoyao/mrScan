#' @title Conduct bidirection MR between traits
#' @param ex_dat1 Dataframe of instruments for the main exposure or outcome. Output from TwoSampleMR::extract_instruments()
#' @param ex_dat2 Dataframe of instruments for a candidate confounder trait. Output from TwoSampleMR::extract_instruments()
#' @param method MR method. See TwoSampleMR::mr_method_list() for details. Default = "mr_raps"
#' @param over.dispersion Parameter for mr_raps method. Default = TRUE
#' @param loss.function Parameter for mr_raps method. Default = "tukey"
#' @param min_instruments minimum number of instruments for candidate traits. Default = 3
#' @returns A list contain bidirection estimates and traits correlation
#'
#' @import TwoSampleMR
#' @import dplyr
#' @export
bidirection_mr <- function(ex_dat1,ex_dat2,method = "mr_raps",
                           over.dispersion = TRUE,loss.function = "tukey",
                           min_instruments=3){
  ID1 <- unique(ex_dat1$id.exposure)
  ID2 <- unique(ex_dat2$id.exposure)
  if(is.null(ex_dat2)){
    return(NULL)
  }else if(!is.null(ex_dat2) & nrow(ex_dat2) < min_instruments){
    return(NULL)
  }else{
    out_dat1 <- extract_outcome_data(snps = ex_dat1$SNP,outcomes = ID2)
    out_dat2 <- extract_outcome_data(snps = ex_dat2$SNP,outcomes = ID1)
    dat_1_2 <- harmonise_data(ex_dat1, out_dat1)
    dat_2_1 <- harmonise_data(ex_dat2, out_dat2)
    m_1_2 <- mr(dat_1_2, method_list = method,
                parameters = list(over.dispersion=over.dispersion,
                                  loss.function=loss.function,
                                  shrinkage=FALSE))
    m_2_1 <- mr(dat_2_1, method_list = method,
                parameters = list(over.dispersion=over.dispersion,
                                  loss.function=loss.function,
                                  shrinkage=FALSE))
    X_1_2 <- dat_1_2 %>%
      rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      select(SNP, beta1, beta2)
    X_2_1 <- dat_2_1 %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    cor_vals<- data.frame(id1 = ID1,id2 = ID2, cor = with(X, cor(beta1, beta2)))
    return(list(mr12 = m_1_2, mr21 = m_2_1, cor = cor_vals))
  }
}
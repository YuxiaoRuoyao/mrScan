#' @title Use ESMR to do MVMR analysis by locally data
#' @param dat A data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for z-values, `trait_ID.ss` for sample sizes.
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.1
#' @returns A dataframe of result summary
#'
#' @import dplyr
#' @import esmr
#' @importFrom purrr map_dfc
#' @export
MVMR_ESMR <- function(dat,R_matrix,pval_threshold=5e-8,effect_size_cutoff=0.1){
  snp <- dat$snp
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
  #beta_hat <- data.frame(beta_hat[, o],check.names = F)
  #se <- data.frame(se[, o],check.names = F)
  z.norm <- z.norm[,o]
  se.norm <- se.norm[,o]
  #i <- ncol(beta_hat)
  i <- ncol(z.norm)
  res.summary <- data.frame(exposure = colnames(z.norm)[-1],
                            b = NA, se = NA, pvalue = NA,
                            method = paste0("ESMR_",pval_threshold))
  filtered_idx <- which(rowSums(abs(data.frame(z.norm[,-1])) < effect_size_cutoff) == ncol(z.norm)-1)
  tryCatch({
    fit <- esmr(beta_hat_Y = z.norm[filtered_idx,1],
                se_Y = se.norm[filtered_idx,1],
                beta_hat_X = z.norm[filtered_idx,2:i],
                se_X = se.norm[filtered_idx, 2:i],
                R = R_matrix,
                pval_thresh = pval_threshold)
    res.summary <- data.frame(exposure = colnames(z.norm)[-1],
                              b = fit$beta$beta_m, se = fit$beta$beta_s) %>%
      mutate(pvalue = 2*pnorm(-abs(b/se)),
             method = paste0("ESMR_",pval_threshold))
  }, error = function(e){
    message("Error in MVMR_ESMR: ", e$message)
  })
  return(res.summary)
}

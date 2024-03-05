#' @title Use GRAPPLE to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats, `trait_ID.se` for standard errors, `trait_ID.p` for pvalues.
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 1e-5
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @returns A dataframe of result summary
#'
#' @import GRAPPLE
#' @import dplyr
#' @import stringr
#' @importFrom purrr map_dfc
#' @export
MVMR_GRAPPLE <- function(dat,R_matrix,pval_threshold = 1e-5,type){
  if(type == "local"){
    beta_hat <- dat %>% select(ends_with(".beta"))
    se <- dat %>% select(ends_with(".se"))
    p <- dat %>% select(ends_with(".p"))
    nms <- stringr::str_replace(names(beta_hat), ".beta", "")
    names(beta_hat)<-names(se)<-names(p)<-nms
    o <- match(colnames(R_matrix), nms)
    beta_hat <- data.frame(beta_hat[, o],check.names = F)
    se <- data.frame(se[, o],check.names = F)
    i <- ncol(beta_hat)
    grapple_dat <- data.frame(cbind(beta_hat, se))
    names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(i-1)),
                            "se_out", paste0("se_exp", 1:(i-1)))
    grapple_dat$selection_pvals <- apply(p[,-1, drop = F],1, min)
    res_and_warning <- withCallingHandlers({
      res <- grappleRobustEst(
        data = grapple_dat,
        plot.it = FALSE,
        p.thres = pval_threshold,
        cor.mat = R_matrix
      )
      res
    }, warning = function(w) {
      if (grepl("Did not converge", w$message)) {
        notConverge <<- TRUE
      } else {
        notConverge <<- FALSE
      }
      invokeRestart("muffleWarning")
    })
    if (!exists("notConverge")) {
      notConverge <- FALSE
    }
    if(i > 2){
      res.summary <- data.frame(exposure=colnames(beta_hat)[-1],
                                b=res_and_warning$beta.hat,
                                se=sqrt(diag(res_and_warning$beta.var)),
                                pvalue=res_and_warning$beta.p.value,
                                method = paste0("GRAPPLE_",pval_threshold),
                                converge = !notConverge)
    }else{
      res.summary <- data.frame(exposure=colnames(beta_hat)[-1],
                                b=res_and_warning$beta.hat,
                                se=sqrt(res_and_warning$beta.var),
                                pvalue=res_and_warning$beta.p.value,
                                method = paste0("GRAPPLE_",pval_threshold),
                                converge = !notConverge)
    }
  }
  if(type == "IEU"){
    id.exposure <- colnames(dat$exposure_beta)
    id.outcome <- dat$outname$id.outcome
    ss.exposure <- gwasinfo(id.exposure)$sample_size
    ss.outcome <- gwasinfo(id.outcome)$sample_size
    z.exposure<- dat$exposure_beta/dat$exposure_se
    z.norm.exposure <- sweep(z.exposure,2,sqrt(ss.exposure),`/`)
    names(ss.exposure) <- id.exposure
    se.norm.exposure <- map_dfc(ss.exposure, ~ rep(1/sqrt(.x),length.out = nrow(dat$exposure_se))) %>%
      as.matrix()
    grapple_dat<-cbind(z.norm.exposure,se.norm.exposure,
                       (dat$outcome_beta/dat$outcome_se)/sqrt(ss.outcome),
                       rep(1/sqrt(ss.outcome),length(dat$outcome_se)))
    i <- length(id.exposure)
    colnames(grapple_dat)<-c(paste0("gamma_exp",seq(1,i)),
                              paste0("se_exp",seq(1,i)),"gamma_out1","se_out1")
    res_and_warning <- withCallingHandlers({
      res <- grappleRobustEst(
        data = grapple_dat,
        plot.it = FALSE
      )
      res
    }, warning = function(w) {
      if (grepl("Did not converge", w$message)) {
        notConverge <<- TRUE
      } else {
        notConverge <<- FALSE
      }
      invokeRestart("muffleWarning")
    })
    if (!exists("notConverge")) {
      notConverge <- FALSE
    }
    res_summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                              b = res_and_warning$beta.hat,
                              se = sqrt(diag(res_and_warning$beta.var)),
                              pvalue = res_and_warning$beta.p.value,
                              method = "MVMR_GRAPPLE",
                              converge = !notConverge)
    res_summary[which(res_summary$se > 1),"pvalue"] <- 1
  }
  return(res.summary)
}

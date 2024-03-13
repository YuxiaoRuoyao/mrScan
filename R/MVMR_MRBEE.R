#' @title Use MRBEE to do MVMR analysis by locally data
#' @param dat If type == "local", it's a data frame of combined GWAS summary data after LD pruning.
#' The required columns include `SNP` for rsID, `trait_ID.beta` for beta hats,
#' `trait_ID.se` for standard errors, `trait_ID.z` for zvalues, `trait_ID.p` for pvalues,
#' `trait_ID.ss` for sample sizes
#' If type == "IEU", it should be the output from TwoSampleMR::mv_harmonise_data().
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param pleio_threshold pvalue threshold in pleiotropy detection. Default = 0
#' @param type Input data type. It could be either "local" for local GWAS summary data after LD pruning or
#' "IEU" for output from TwoSampleMR::mv_harmonise_data().
#' @param ss.exposure A vector of sample size for exposures. You can provide it when type = "IEU".
#' The order of it should be the same with beta hat matrix and se matrix. Default = NULL
#' @param ss.outcome A numeric number of sample size for the outcome. You can provide it when type = "IEU". Dafault = NULL
#' @param effect_size_cutoff Standardized effect size threshold. Default = 0.05
#' @returns A dataframe of result summary
#'
#' @import MRBEE
#' @import dplyr
#' @import ieugwasr
#' @importFrom purrr map_dfc
#' @export
MVMR_MRBEE <- function(dat,R_matrix,pval_threshold = 5e-8,pleio_threshold = 0,type,
                       ss.exposure = NULL, ss.outcome = NULL, effect_size_cutoff=0.05){
  if(type == "local"){
    z <- dat %>% select(ends_with(".z"))
    p <- dat %>% select(ends_with(".p"))
    ss <- dat %>% select(ends_with(".ss"))
    nms <- stringr::str_replace(names(z), ".z", "")
    #names(p)<-names(z)<-nms
    names(z)<-names(p)<-names(ss)<-nms
    N <- apply(ss, 2, median, na.rm = TRUE)
    z.norm <- sweep(z,2,sqrt(N),`/`) %>% data.frame(check.names = F)
    se.norm <- purrr::map_dfc(ss, ~ rep(1/sqrt(.x),length.out = nrow(z))) %>%
      data.frame(check.names = F)
    o <- match(colnames(R_matrix), nms)
    #z <- data.frame(z[, o],check.names = F)
    z.norm <- z.norm[,o]
    se.norm <- se.norm[,o]
    #i <- ncol(z)
    i <- ncol(z.norm)
    pmin <- apply(p[,-1, drop = F], 1, min)
    ix <- which(pmin < pval_threshold)
    # Make the last one be outcome for R matrix
    R_matrix <- R_matrix[c(nms[-1],nms[1]),c(nms[-1],nms[1])]
    if(i>2){
      # fit <- MRBEE.IMRP(by=z[ix,1],bX=as.matrix(z[ix,-1]),
      #                   byse=rep(1,length(ix)),
      #                   bXse=matrix(1,length(ix),i-1),
      #                   Rxy=R_matrix,
      #                   pv.thres = pleio_threshold, var.est = "variance")
      fit <- MRBEE.IMRP(by=z.norm[ix,1],bX=as.matrix(z.norm[ix,-1]),
                        byse=se.norm[ix,1],
                        bXse=as.matrix(se.norm[ix,-1]),
                        Rxy=R_matrix,
                        pv.thres = pleio_threshold, var.est = "variance")
      res.summary <- data.frame(exposure = names(fit$theta),
                                b = fit$theta,
                                se = sqrt(diag(fit$covtheta)))
    }else{
      # fit <- MRBEE.IMRP.UV(by = z[ix,1],bx = z[ix,-1],
      #                      byse = rep(1,length(ix)),
      #                      bxse = rep(1,length(ix)),
      #                      Rxy=R_matrix,
      #                      pv.thres = pleio_threshold, var.est="variance")
      fit <- MRBEE.IMRP.UV(by = z.norm[ix,1],bx = z.norm[ix,-1],
                           byse = se.norm[ix,1],
                           bxse = se.norm[ix,-1],
                           Rxy=R_matrix,
                           pv.thres = pleio_threshold, var.est="variance")
      res.summary <- data.frame(exposure = nms[-1],
                                b = fit$theta,
                                se = sqrt(fit$vartheta))
    }
    res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
    res.summary$method <- paste0("MRBEE_",pval_threshold,"_pleio_",pleio_threshold)
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
    filtered_idx <- which(rowSums(abs(z.norm.exposure) < effect_size_cutoff) == ncol(z.norm.exposure))

    fit <- MRBEE.IMRP(by = z.norm.outcome[filtered_idx],
                      bX = z.norm.exposure[filtered_idx, ],
                      byse = se.norm.outcome[filtered_idx],
                      bXse = se.norm.exposure[filtered_idx, ],
                      Rxy = diag(nrow = length(id.exposure)+1),
                      pv.thres = pleio_threshold, var.est = "variance")
    res.summary <- data.frame(id.exposure = id.exposure, id.outcome = id.outcome,
                              b = fit$theta,se = sqrt(diag(fit$covtheta))) %>%
      mutate(pvalue = 2*pnorm(-abs(b/se)), method = "MVMR_MRBEE")
  }
  res.summary[which(res.summary$se > 1),"pvalue"] <- 1
  rownames(res.summary) <- NULL
  return(res.summary)
}

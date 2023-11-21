#' @import dplyr
#' @importFrom dplyr select filter rename
#' @import TwoSampleMR
#' @import mr.raps
#' @import plyr
#' @import purrr
make_metadata <- function(select_id,id,trait,sex,consortium=NA,nsnp,unit=NA,author=NA,
                          note=NA,sample_size,pmid=NA,population,year=NA,category=NA,
                          doi=NA,sd=NA,ncase=NA,ncontrol=NA){
  df <- data.frame(id=id,trait=trait,sex=sex,consortium=consortium,
             nsnp=nsnp,unit=unit,author=author,note=note,
             sample_size=sample_size,pmid=pmid,population=population,year=year,
             category=category,doi=doi,sd=sd,ncase=ncase,ncontrol=ncontrol)
  # select_id is result from extract_traits$id.list
  df <- df %>% filter(id %in% select_id) %>% mutate(status = "Initial List")
  return(df)
}
mr_new <- function (dat, parameters = default_parameters(), method_list = subset(mr_method_list(),
                                                                                 use_by_default)$obj)
{
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),
                        function(x1) {
                          x <- subset(x1, mr_keep)
                          if (nrow(x) == 0) {
                            message("No SNPs available for MR analysis of '",
                                    x1$id.exposure[1], "' on '", x1$id.outcome[1],
                                    "'")
                            return(NULL)
                          }
                          else {
                            message("Analysing '", x1$id.exposure[1], "' on '",
                                    x1$id.outcome[1], "'")
                          }
                          res <- lapply(c("mr_raps"), function(meth) {
                            get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure,
                                      x$se.outcome, parameters)
                          })
                          methl <- mr_method_list()
                          mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1],
                                               method = methl$name[match(method_list, methl$obj)],
                                               nsnp = sapply(res, function(x) x$nsnp), b = sapply(res,
                                                                                                  function(x) x$b), se = sapply(res, function(x) x$se),
                                               pval = sapply(res, function(x) x$pval))
                          mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) &
                                                       is.na(pval)))
                          return(mr_tab)
                        })
  return(mr_tab)
}
mr_raps<-function (b_exp, b_out, se_exp, se_out,parameters = default_parameters()) {
  data <- data.frame(beta.exposure = b_exp, beta.outcome = b_out,
                     se.exposure = se_exp, se.outcome = se_out)
  out <- suppressWarnings(mr.raps::mr.raps(b_exp,b_out, se_exp, se_out,
                                           diagnosis = FALSE,
                                           over.dispersion = parameters$over.dispersion,
                                           loss.function = parameters$loss.function))
  list(b = out$beta.hat, se = out$beta.se,
       pval = stats::pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
}
mr_pairs_raps<-function (ids1, ids2, inst_pval = 5e-08,
                         over.dispersion=TRUE,loss.function = "tukey"){
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids1, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP,
                                               outcomes = ids2)
  dat_1_2 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_1_2 <- mr_new(dat_1_2, method_list = c("mr_raps"),
              parameters = list(over.dispersion=over.dispersion,
                                loss.function=loss.function))
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_2_1 <- mr_new(dat_2_1, method_list = c("mr_raps"),
                  parameters = list(over.dispersion=over.dispersion,
                                    loss.function=loss.function))
  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
    filter(id1 != id2)
  cor_vals$cor <- purrr::map2(cor_vals$id1, cor_vals$id2, function(x, y) {
    X_1_2 <- dplyr::filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
      dplyr::rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      dplyr::select(SNP, beta1, beta2)
    X_2_1 <- dplyr::filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      dplyr::rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      dplyr::select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    with(X, cor(beta1, beta2))
  }) %>% unlist()
  return(list(mr12 = m_1_2, mr21 = m_2_1, cor = cor_vals))
}
calculate_cor <- function (ids1, ids2, inst_pval = 5e-08){
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids1, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP,
                                               outcomes = ids2)
  dat_1_2 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
    filter(id1 != id2)
  cor_vals$cor <- purrr::map2(cor_vals$id1, cor_vals$id2, function(x, y) {
    X_1_2 <- filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
      rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      select(SNP, beta1, beta2)
    X_2_1 <- filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    with(X, cor(beta1, beta2))
  }) %>% unlist()
  return(cor_vals)
}
run_grapple <- function(beta.exposure,beta.outcome,se.exposure,se.outcome){
  grapple.data<-data.frame(cbind(beta.outcome, beta.exposure,
                                 se.outcome, se.exposure))
  i <- ncol(beta.exposure)
  names(grapple.data)<-c("gamma_out", paste0("gamma_exp", 1:i),
                         "se_out", paste0("se_exp", 1:i))
  res <- GRAPPLE::grappleRobustEst(data = grapple.data,
                          plot.it =FALSE,
                          niter = 100000)
  res.summary <- data.frame(exposure=colnames(beta.exposure),
                            b=res$beta.hat,
                            se=sqrt(diag(res$beta.var)),
                            pvalue=res$beta.p.value,
                            method = "GRAPPLE")
  return(res.summary)
}
run_MRBEE <- function(beta.exposure,beta.outcome,se.exposure,se.outcome,
                      pleio_p_thresh = 0, Rcor = NULL){
  beta_hat <- data.frame(cbind(beta.outcome,beta.exposure))
  se <- data.frame(cbind(se.outcome,se.exposure))
  if(is.null(Rcor)){
    Rcor <- diag(1,nrow = ncol(beta_hat))
    print("Assume independence between traits!")
  }
  bT <- list(R = Rcor, Ncor = Inf,
             EstHarm = beta_hat,
             SEHarm =  se)
  pD <- MRBEE::prepData(bT,verbose =FALSE)
  fit <- MRBEE::MRBEE.IMRP(pD, PleioPThreshold = pleio_p_thresh)
  res.summary <- data.frame(exposure = colnames(beta.exposure),
                            b = fit$CausalEstimates[-1],
                            se = sqrt(diag(fit$VCovCausalEstimates))[-1])
  res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
  res.summary$method <- "MRBEE"
  return(res.summary)
}

#' @import dplyr
#' @importFrom dplyr select filter rename group_by summarize
#' @import TwoSampleMR
#' @import ieugwasr
#' @import mr.raps
#' @import plyr
#' @import purrr
retrieve_traits <- function (id_x, pval_x = 5e-8, pval_z = 1e-5,
                             pop = "EUR", batch = c("ieu-a", "ieu-b","ukb-b"),
                             r2 = 0.001, kb = 10000,
                             access_token = check_access_token(), min_snps = 5) {
  top_hits <- tophits(id = id_x, pval = pval_x, r2 = r2, kb = kb,
                      pop = pop, access_token = access_token)
  cat("Retrieved", nrow(top_hits), "instruments for", id_x,
      "\n")
  phe <- ieugwasr::phewas(variants = top_hits$rsid, pval = pval_z, batch = batch,
                access_token = access_token)
  cat("Retrieved", nrow(phe), "associations with", length(unique(phe$id)),
      "traits", "\n")
  x <- phe %>% dplyr::group_by(id) %>% dplyr::summarize(n = length(unique(rsid)))
  cat(sum(x$n >= min_snps), "traits have at least", min_snps,
      "shared variants with", id_x, "\n")
  ids <- x$id[x$n >= min_snps]
  phe <- dplyr::filter(phe, id %in% ids)
  ret <- list(topx = top_hits, phe = phe)
  return(ret)
}
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
mr_pairs<-function (ids1, ids2, inst_pval = 5e-08,method_list = c("mr_raps"),
                    over.dispersion=TRUE,loss.function = "tukey",shrinkage=FALSE){
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids1, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP,
                                               outcomes = ids2)
  dat_1_2 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_1_2 <- mr(dat_1_2, method_list = method_list,
              parameters = list(over.dispersion=over.dispersion,
                                loss.function=loss.function,
                                shrinkage=shrinkage))
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_2_1 <- mr(dat_2_1, method_list = method_list,
              parameters = list(over.dispersion=over.dispersion,
                                loss.function=loss.function,
                                shrinkage=shrinkage))
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
  # This step redundant
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
    filter(id1 != id2)
  cor_vals$cor <- purrr::map2(cor_vals$id1, cor_vals$id2, function(x, y) {
    X_1_2 <- dplyr::filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
      rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      select(SNP, beta1, beta2)
    X_2_1 <- dplyr::filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- dplyr::bind_rows(X_1_2, X_2_1)
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
download_gwas <- function(id_list,position=NULL){
  f1 <- paste0("wget https://gwas.mrcieu.ac.uk/files/",id_list,"/",id_list,".vcf.gz")
  f2 <- paste0("wget https://gwas.mrcieu.ac.uk/files/",id_list,"/",id_list,".vcf.gz.tbi")
  f <- data.frame(c(f1,f2))
  write.table(f,file = paste0(position,"download.sh"),row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}

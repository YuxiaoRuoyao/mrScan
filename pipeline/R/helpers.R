library(plyr)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(mr.raps)
library(purrr)
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
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
    filter(id1 != id2)
  cor_vals$cor <- purrr::map2(cor_vals$id1, cor_vals$id2, function(x, y) {
    X_1_2 <- dplyr::filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
      dplyr::rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
      dplyr::select(SNP, beta1, beta2)
    X_2_1 <- dplyr::filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      dplyr::rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      dplyr::select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
    X <- dplyr::bind_rows(X_1_2, X_2_1)
    with(X, cor(beta1, beta2))
  }) %>% unlist()
  return(cor_vals)
}
# calculate_cor <- function (df_pairs, inst_pval = 5e-08){
#   df_pairs$cor <- purrr::map2(df_pairs$X1, df_pairs$X2, function(x, y){
#     ex_dat <- extract_instruments(outcomes = x, p1 = inst_pval)
#     out_dat <- extract_outcome_data(snps = ex_dat$SNP,outcomes = y)
#     dat_1_2 <- harmonise_data(ex_dat, out_dat)
#     ex_dat <- extract_instruments(outcomes = y, p1 = inst_pval)
#     out_dat <- extract_outcome_data(snps = ex_dat$SNP, outcomes = x)
#     dat_2_1 <- harmonise_data(ex_dat, out_dat)
#     X_1_2 <- dplyr::filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
#       rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
#       select(SNP, beta1, beta2)
#     X_2_1 <- dplyr::filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
#       rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
#       select(SNP, beta1, beta2) %>% filter(!SNP %in% X_1_2$SNP)
#     X <- dplyr::bind_rows(X_1_2, X_2_1)
#     with(X, cor(beta1, beta2))
#   }) %>% unlist()
#   return(df_pairs)
# }
calculate_cor_pairwise <- function(id.list,inst_pval = 5e-8){
  my_inst <- purrr::map(id.list, function(x){
    df_inst <- extract_instruments(outcomes = x, p1 = inst_pval)
    print(paste0("Finish ",x," !"))
    return(df_inst$SNP)
  })
  all_inst <- unique(unlist(my_inst))
  chunk_inst <- split(all_inst, ceiling(seq_along(all_inst)/20))
  beta_matrix <- purrr::map_dfc(id.list,function(t){
    f <- purrr::map_dfr(chunk_inst, function(x) ieugwasr::associations(x,t) %>%
             distinct(rsid, .keep_all = TRUE)) # Ask this: how to deal with duplicated items?
    new_f <- data.frame(SNP = all_inst) %>%
      left_join(f[,c("rsid","beta")],by=c("SNP" = "rsid")) %>%
      select(beta)
    colnames(new_f) <- t
    print(paste0("Finish ",t," !"))
    return(new_f)
  })
  beta_matrix <- cbind(data.frame(SNP=all_inst),beta_matrix)
  names(my_inst) <- id.list
  res <- data.frame(t(combn(id.list,2)))
  res$cor <- map2(res$X1, res$X2, function(i, j) {
    sub_inst <- my_inst[c(i, j)]
    inst <- unique(unlist(sub_inst))
    sub_matrix <- beta_matrix %>% filter(SNP %in% inst)
    cor(sub_matrix[[i]],sub_matrix[[j]],use = "complete.obs")
  }) %>% unlist()
  res2 <- data.frame(X1 = res$X2, X2 = res$X1, cor = res$cor)
  res_all <- rbind(res,res2)
  df_matrix <- as.data.frame.matrix(xtabs(cor ~ ., res_all))
  df_matrix <- abs(df_matrix)
  return(list(df_pair = res, R_matrix = df_matrix))
}
run_grapple <- function(beta.exposure,beta.outcome,se.exposure,se.outcome,R=NULL){
  grapple.data<-data.frame(cbind(beta.outcome, beta.exposure,
                                 se.outcome, se.exposure))
  i <- ncol(beta.exposure)
  names(grapple.data)<-c("gamma_out", paste0("gamma_exp", 1:i),
                         "se_out", paste0("se_exp", 1:i))
  res <- GRAPPLE::grappleRobustEst(data = grapple.data,
                          plot.it =FALSE,cor.mat = R,
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
format_ieu_chrom <- function(file, chrom){
  dat <- query_chrompos_file(paste0(chrom, ":1-536870911"), file) %>%
    vcf_to_tibble() %>%
    mutate(p_value = 10^{-1*LP})
  dat <- sumstatFactors::gwas_format(dat, "ID", "ES", "SE", "ALT",
                     "REF", "seqnames", "start",
                     p_value = "p_value",
                     sample_size = "SS",
                     compute_pval = TRUE)
  return(dat)
}
format_flat_chrom <- function(file, chrom,
                              snp_name,
                              pos_name,
                              chrom_name,
                              A1_name, A2_name,
                              beta_hat_name,
                              se_name,
                              p_value_name,
                              af_name,
                              sample_size_name,
                              effect_is_or
){
  if(!p_value_name %in% c("", "NA", NA)){
    pstring <- paste0(", `", p_value_name, "`='d'")
  }else{
    pstring <- ""
    p_value_name <- NA
  }
  if(!sample_size_name %in% c("", "NA", NA)){
    sstring <- paste0(", `", sample_size_name, "`='d'")
  }else{
    sstring <- ""
    sample_size_name <- NA
  }
  if(!pos_name %in% c("", "NA", NA)){
    posstring <- paste0(", `", pos_name, "`='d'")
  }else{
    posstring <- ""
    pos_name <- NA
  }
  if(!af_name %in% c("", "NA", NA)){
    afstring <- paste0(", `", af_name, "`='d'")
  }else{
    afstring <- ""
    af_name <- NA
  }
  col_string <- paste0("cols_only(`", snp_name, "`='c', `",
                       A1_name , "`='c', `", A2_name, "`='c', `",
                       beta_hat_name , "`='d', `", se_name, "`='d', `",
                       chrom_name, "`='c' ", posstring,
                       pstring,  sstring, afstring, ")")


  if(str_ends(file, "gz") ){
    h <- read_table(pipe(paste0("gzip -cd ", file, " | head -2")))
    n <- which(names(h) == chrom_name)
    awk_cmd <- paste0("gzip -cd ", file, " | awk '{if ($", n, " == \"", chrom, "\") print $0}' - ")
  }else{
    h <- read_table(pipe(paste0("head -2 ", file)))
    n <- which(names(h) == chrom_name)
    awk_cmd <- paste0("awk '{if ($", n, " == \"", chrom, "\") print $0}' ", file)
  }
  X <- read_table(pipe(awk_cmd), col_types = eval(parse(text = col_string)), col_names = names(h))

  if(!is.na(af_name)){
    ix <- which(X[[af_name]] > af_thresh & X[[af_name]] < (1-af_thresh))
    X <- X[ix,]
  }

  if(effect_is_or){
    X$beta <- log(X[[beta_hat_name]])
    beta_hat <- "beta"
  }
  dat <- gwas_format(X, snp_name, beta_hat_name, se_name, A1_name,
                     A2_name, chrom_name, pos_name,
                     p_value = p_value_name,
                     sample_size = sample_size_name,
                     compute_pval = TRUE)
  return(dat)
}

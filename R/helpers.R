#' @import dplyr
#' @import TwoSampleMR
#' @import ieugwasr
#' @import mr.raps
#' @import stringr
#' @import sumstatFactors
#' @import gwasvcf
#' @import rlang
#' @import readr
#' @importFrom purrr map_dfr map2 map
#' @export
get_inst_IEU <- function(id_x,pval_x = 5e-8,r2 = 0.001, kb = 10000, pop = "EUR",
                         access_token = ieugwasr::check_access_token()){
  top_hits <- tophits(id = id_x, pval = pval_x, r2 = r2, kb = kb,
                      pop = pop, access_token = access_token)
  cat("Retrieved", nrow(top_hits), "instruments for", id_x,
      "\n")
  return(top_hits)
}
# EBI function is not correct now.
get_inst_EBI <- function(id_x, pval_x = 5e-8){
  df_associations <- gwasrapidd::get_associations(study_id = id_x,warnings = FALSE)
  tbl01 <- dplyr::select(df_associations@risk_alleles, association_id, variant_id, risk_allele)
  tbl02 <- dplyr::select(df_associations@associations, association_id, pvalue, beta_number, beta_unit, beta_direction)
  df_variants <- dplyr::left_join(tbl01, tbl02, by = 'association_id') %>%
    tidyr::drop_na() %>%
    dplyr::arrange(variant_id, risk_allele) %>%
    dplyr::filter(pvalue < pval_x) %>%
    dplyr::rename(rsid = variant_id)
  cat("Retrieved", nrow(df_variants), "instruments for", id_x,
      "\n")
  return(df_variants)
}
#' @export
get_inst_local <- function(file_path,id_x,r2 = 0.001, kb = 10000,ref_path,pval_x=5e-8,
                           snp_name = NA, pos_name = NA, chrom_name = NA, A1_name = NA,
                           A2_name = NA, beta_hat_name = NA, se_name = NA, p_value_name = NA,
                           af_name = NA, sample_size_name = NA, effect_is_or = NA){
  fulldat <- purrr::map_dfr(seq(1:22), function(c){
    if(str_ends(file_path, "vcf.gz") | str_ends(file_path, "vcf.bgz")){
      dat <- format_ieu_chrom(file_path, c)
    }else if(str_ends(file_path, ".h.tsv.gz")){
      dat <- format_flat_chrom(file_path, c,
                               snp_name = "hm_rsid",
                               pos_name = "hm_pos",
                               chrom_name = "hm_chrom",
                               A1_name = "hm_effect_allele",
                               A2_name = "hm_other_allele",
                               beta_hat_name = "hm_beta",
                               se_name = "standard_error",
                               p_value_name = "p_value",
                               af_name = "hm_effect_allele_frequency",
                               sample_size_name = NA,
                               effect_is_or = FALSE)
    }else{
      if(is.na(beta_hat_name)|is.na(snp_name)|is.na(chrom_name)|is.na(se_name)|is.na(A1_name)|is.na(A2_name)){
        stop("SNP, Beta_hat, SE, A1 and A2 column names are required.\n")
      }
      dat <- format_flat_chrom(file_path, c,
                               snp_name = snp_name,
                               pos_name = pos_name,
                               chrom_name = chrom_name,
                               A1_name = A1_name,
                               A2_name = A2_name,
                               beta_hat_name = beta_hat_name,
                               se_name = se_name,
                               p_value_name = p_value_name,
                               af_name = af_name,
                               sample_size_name = sample_size_name,
                               effect_is_or = effect_is_or)
    }
    pos_name <- as_name(paste0(id_x, ".pos"))
    beta_name <- as_name(paste0(id_x, ".beta"))
    se_name <- as_name(paste0(id_x, ".se"))
    p_name <- as_name(paste0(id_x, ".p"))
    z_name <- as_name(paste0(id_x, ".z"))
    dat <- dat %>% dplyr::mutate(Z = beta_hat/se) %>%
      dplyr::rename(REF = A2, ALT = A1) %>%
      dplyr::select(chrom, snp, REF, ALT,
                    !!pos_name := pos,
                    !!beta_name := beta_hat,
                    !!se_name := se,
                    !!p_name := p_value,
                    !!z_name := Z)
    dup_snps <- dat$snp[duplicated(dat$snp)]
    if(length(dup_snps) > 0){
      dat <- filter(dat, !snp %in% dup_snps)
    }
    ld_prune_plink(X = dat,r2_thresh = r2,clump_kb = kb,ref_path = ref_path)
  })
  p <- fulldat %>% select(ends_with(".p"))
  ix <- which(p < pval_x)
  df_inst <- fulldat[ix,] %>%
    dplyr::rename_with(~ c("chr","rsid","position"), all_of(c("chrom","snp",paste0(id_x, ".pos"))))
  cat("Retrieved", nrow(df_inst), "instruments for", id_x,"\n")
  return(df_inst)
}
#' @export
get_exposure_inst <- function(id_x,type,file_path = NA,pval_x = 5e-8,r2 = 0.001,
                              kb = 10000, pop = "EUR",
                              access_token = ieugwasr::check_access_token(),
                              ref_path = NA){
  if(type == "IEU"){
    df_inst <- get_inst_IEU(id_x = id_x,pval_x = pval_x,r2 = r2,kb = kb,pop = pop,
                 access_token = access_token)
  }else if(type == "local"){
    if(is.na(file_path)){
      print("Please provide local file path!")
    }else if(is.na(ref_path)){
      print("Please provide LD reference file path!")
    }
    df_inst <- get_inst_local(file_path = file_path,id_x = id_x,r2 = r2,kb = kb,
                              pval_x = pval_x,ref_path = ref_path)
  }
  return(df_inst)
}
#' @export
get_association_IEU <- function(df_inst,pval_z = 1e-5,batch = c("ieu-a", "ieu-b","ukb-b"),
                                access_token = ieugwasr::check_access_token()){
  counter <- 0
  progress_message <- function() {
    counter <<- counter + 1
    if (counter %% 50 == 0) {
      cat("Processed", counter, "SNPs...\n")
    }
  }
  results <- purrr::map(df_inst$rsid, function(rsid) {
    progress_message()
    ieugwasr::phewas(variants = rsid, pval = pval_z, batch = batch,
                     access_token = access_token)
  })
  phe <- dplyr::bind_rows(results)
  cat("Retrieved", nrow(phe), "associations with", length(unique(phe$id)),
      "traits", "\n")
  return(phe)
}

# get_association_IEU <- function(df_inst,pval_z = 1e-5,batch = c("ieu-a", "ieu-b","ukb-b"),
#                                 access_token = ieugwasr::check_access_token()){
#   batch1 <- batch[batch %in% c("ieu-a", "ieu-b","ukb-b")]
#   batch2 <- batch[!batch %in% c("ieu-a", "ieu-b","ukb-b")]
#   phe1 <- ieugwasr::phewas(variants = df_inst$rsid, pval = pval_z, batch = batch1,
#                            access_token = access_token)
#   if(length(batch2) != 0){
#     phe2 <- ieugwasr::phewas(variants = df_inst$rsid, pval = pval_z, batch = batch2,
#                              access_token = access_token)
#     phe <- rbind(phe1,phe2)
#   }else{
#     phe <- phe1
#   }
#   cat("Retrieved", nrow(phe), "associations with", length(unique(phe$id)),
#       "traits", "\n")
#   return(phe)
# }
#' @export
get_association_local <- function(file_path,trait_id,df_inst,pval_z,snp_name=NA,
                                  beta_hat_name = NA, se_name = NA,
                                  p_value_name = NA){
  tmp_file <- tempfile()
  writeLines(df_inst$rsid, tmp_file)
  if(str_ends(file_path, "vcf.gz") | str_ends(file_path, "vcf.bgz")){
    dat_filter <- query_chrompos_file(chrompos = paste0(df_inst$chr,":",df_inst$position,"-",df_inst$position),
                                      vcffile = file_path) %>%
      vcf_to_tibble() %>%
      mutate(p_value = 10^{-1*LP}) %>%
      filter(p_value < pval_z)
  }else if(str_ends(file_path, ".h.tsv.gz")){
    h <- readr::read_table(pipe(paste0("gzip -cd ", file_path, " | head -2")))
    n <- which(names(h) == "hm_rsid")
    awk_cmd <- paste0("gzip -cd ", file_path, " | awk 'NR==FNR { A[$1]=1 ; next } $",
                      n ," in A' ", tmp_file, " -")
    dat_filter <- read_table(pipe(awk_cmd), col_names = names(h)) %>%
      dplyr::rename(rsid = hm_rsid) %>%
      mutate(p_value = 2 * pnorm(-abs(hm_beta/standard_error))) %>%
      filter(p_value < pval_z)
  }else{
    if(is.na(snp_name) | is.na(beta_hat_name) | is.na(se_name)){
      stop("Please input the column names of SNP, Beta hat and SE!")
    }
    if(str_ends(file_path, ".gz")){
      h <- readr::read_table(pipe(paste0("gzip -cd ", file_path, " | head -2")))
      n <- which(names(h) == snp_name)
      awk_cmd <- paste0("gzip -cd ", file_path, " | awk 'NR==FNR { A[$1]=1 ; next } $",
                        n ," in A' ", tmp_file, " -")
    }else{
      h <- readr::read_table(pipe(paste0("head -2 ", file_path)))
      n <- which(names(h) == snp_name)
      awk_cmd <- paste0("awk 'NR==FNR { A[$1]=1 ; next } $", n ," in A' ",
                        tmp_file, " ", file_path)
    }
    X <- read_table(pipe(awk_cmd), col_names = names(h))
    if(!is.na(p_value_name)){
      dat_filter <- X %>% dplyr::rename(p_value = p_value_name,rsid = snp_name) %>%
        filter(p_value < pval_z)
    }else{
      dat_filter <- X %>% dplyr::rename(beta_hat = beta_hat_name, se = se_name,rsid = snp_name) %>%
        mutate(p_value = 2 * pnorm(-abs(beta_hat/se))) %>%
        filter(p_value < pval_z)
    }
  }
  dat_filter$id <- trait_id
  unlink(tmp_file)
  return(dat_filter)
}
#' @export
get_association_inst <- function(df_inst,type,pval_z = 1e-5,batch = c("ieu-a", "ieu-b","ukb-b"),
                                 access_token = ieugwasr::check_access_token(),
                                 file_list = NA, trait_list = NA,
                                 snp_name_list=NA,beta_hat_name_list = NA,
                                 se_name_list = NA,p_value_name_list = NA){
  if(type == "IEU"){
    df_association <- get_association_IEU(df_inst = df_inst,pval_z = pval_z,batch = batch,
                                          access_token = access_token)
  }else if(type == "local"){
    df_association <- purrr::map_dfr(seq(length(file_list)),function(i){
      get_association_local(file_path = file_list[i],trait_id = trait_list[i],
                            df_inst = df_inst,pval_z = pval_z,
                            snp_name=snp_name_list[i],
                            beta_hat_name = beta_hat_name_list[i],
                            se_name = se_name_list[i],
                            p_value_name = p_value_name_list[i])
    })
    cat("Retrieved", nrow(df_association), "associations with", length(file_list),
        "traits", "\n")
  }
  return(df_association)
}
#' @export
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
#' @export
download_gwas <- function(id_list,df_harmonise = NULL,data_path = NULL,
                          path_checkpoint = NULL){
  ebi_list <- id_list[grep("GCST",id_list)]
  regular_list <- id_list[!id_list %in% ebi_list]
  if(length(ebi_list)>0 & is.null(df_harmonise)){
    stop("Please provide the path of harmonised_list.txt")
  }else if(length(ebi_list)>0 & !is.null(df_harmonise)){
    GCST_list <- ebi_list %>% strsplit("-") %>% sapply(tail,1)
    df_harmonise$V2 <- df_harmonise$V1 %>% strsplit("-") %>% sapply( "[", 3)
    GCST_file <- df_harmonise %>% filter(str_detect(V2,paste0(GCST_list,collapse = '|'))) %>%
      pull(V1) %>% str_split("/",n=2) %>% sapply(tail,1)
    nf_list <- ebi_list[!GCST_list %in% df_harmonise$V2]
    if(length(nf_list) > 0){
      cat(paste0("Cannot find correct download link of ",nf_list,". Please input manually!"))
    }
    f3 <- paste0("wget -N -P ",data_path," https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",GCST_file)
  }else{
    f3 <- NULL
  }
  f1 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",regular_list,"/",regular_list,".vcf.gz")
  f2 <- paste0("wget -N -P ",data_path," https://gwas.mrcieu.ac.uk/files/",regular_list,"/",regular_list,".vcf.gz.tbi")
  checkpoint <- paste0("echo 'all done!' > ",path_checkpoint)
  f <- data.frame(c(f1,f2,f3,checkpoint))
  return(f)
}
#' @export
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
#' @export
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
#' @export
filter_high_cor_XY <- function(id_list,df_info,res_cor,id_exposure,R2_cutoff){
  trait_cor_X <- res_cor %>% filter(id1 == id_exposure) %>%
    filter(abs(cor) > R2_cutoff) %>%
    pull(id2)
  trait_cor_Y <- res_cor %>% filter(id1 != id_exposure) %>%
    filter(abs(cor) > R2_cutoff) %>%
    pull(id2)
  df_info[df_info$id %in% trait_cor_X,"status"] <- "delete since high cor with X"
  df_info[df_info$id %in% trait_cor_Y,"status"] <- "delete since high cor with Y"
  select_trait <- id_list[!id_list %in% c(trait_cor_X,trait_cor_Y)]
  return(list(id.list = select_trait, trait.info = df_info))
}

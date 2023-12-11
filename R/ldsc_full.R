#' @title Estimate genetic correlation matrix
#' @param beta_dir Directory path of merged data. Default is the current work directory.
#' @param ref_path Path for the LD reference panel.
#' @param out_dir Output data path. Default is in the current work directory.
#' @param prefix Name prefix for the output. Default = NULL
#' @returns Save one dataframe per chromosome with columns for SNP info
#'
#' @import ieugwasr
#' @import dplyr
#' @import purrr
#' @import readr
#' @import bigsnpr
#' @import stringr
#' @export
# Add a pvalue option for calculate R
ldsc_full <- function(l2_dir,beta_dir=NULL,out_dir=NULL,prefix=NULL){
  beta_files <- paste0(beta_dir,prefix,".beta.",seq(1,22),".RDS")
  ld <- purrr::map_dfr(1:22, function(c){
    read_table(paste0(l2_dir, c, ".l2.ldscore.gz"))
  })

  M <- purrr:::map(1:22, function(c){
    read_lines(paste0(l2_dir, c, ".l2.M_5_50"))
  }) %>% unlist() %>% as.numeric() %>% sum()

  dat <- map_dfr(beta_files, function(x){readRDS(x)})
  dat <- filter(dat, snp %in% ld$SNP)

  Z <- dat %>% dplyr::select(ends_with(".z"))
  n <- names(Z)
  n <- stringr::str_replace(n, ".z", "")
  names(Z) <- n
  Z$SNP <- dat$snp

  full_dat <- inner_join(Z, ld)

  h2 <- lapply(n, function(nn){
    ss <- median(dat[[paste0(nn,".ss")]])
    i <- which(!is.na(full_dat[[nn]]))
    snp_ldsc(ld_score = full_dat$L2[i],
             ld_size = M,
             chi2 = (full_dat[[nn]][i])^2,
             sample_size = ss,
             blocks = NULL)
  })

  name_pairs <- expand.grid(t1 = seq_along(n), t2 = seq_along(n)) %>%
    filter(t1 < t2)
  name_pairs$t1 <- n[name_pairs$t1]
  name_pairs$t2 <- n[name_pairs$t2]

  np_res <- lapply(seq(nrow(name_pairs)), function(ii){
    cat(ii, " ")
    nn1 <- name_pairs$t1[ii]
    nn2 <- name_pairs$t2[ii]
    ss1 <- median(dat[[paste0(nn1,".ss")]])
    ss2 <- median(dat[[paste0(nn2,".ss")]])
    i <- which(!is.na(full_dat[[nn1]]) & !is.na(full_dat[[nn2]]))
    snp_ldsc_rg(ld_score = full_dat$L2[i],
                ld_size = M,
                sample_size_1 = ss1,
                sample_size_2 = ss2,
                z1 = full_dat[[nn1]][i],
                z2 = full_dat[[nn2]][i],
                blocks = NULL, h2_se = FALSE)
  })
  name_pairs$cov <- np_res %>% map("int") %>% unlist()
  x <- name_pairs
  names(x) <- c("t2", "t1", "cov")
  name_pairs <- bind_rows(name_pairs, data.frame(t1 = n, t2 = n, cov = h2 %>% map("int") %>% unlist()))
  name_pairs <- bind_rows(name_pairs, x)
  cov_mat <- reshape2::dcast(name_pairs, t1 ~ t2)

  nms <- paste0(as.vector(cov_mat$t1))
  R <- as.matrix(cov_mat[,-1])
  R <- Matrix::nearPD(R,  posd.tol = 1e-3) %>% with(., as.matrix(mat))
  o <- match(n, nms)
  R <- R[o, o]
  Rcor <- cov2cor(R)
  saveRDS(Rcor, file=paste0(out_dir,prefix,".R_est.RDS"))
}

#' @title Estimate genetic correlation matrix by LDSC method
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
#' @import sumstatFactors
#' @export
ldsc_full<-function(l2_dir,beta_dir=NULL,out_dir=NULL,prefix=NULL){
  beta_files <- paste0(beta_dir,prefix,".beta.",seq(1,22),".RDS")
  ld <- purrr::map_dfr(1:22, function(c){
    read_table(paste0(l2_dir, c, ".l2.ldscore.gz"))
  })
  M <- purrr:::map(1:22, function(c){
    read_lines(paste0(l2_dir, c, ".l2.M_5_50"))
  }) %>% unlist() %>% as.numeric() %>% sum()
  X <- map_dfr(beta_files, function(f){
    readRDS(f) %>%
      rename(SNP = snp) %>%
      inner_join(., ld)})
  Z_hat <- X %>%
    select(ends_with(".z")) %>%
    as.matrix()
  nmsz <- str_replace(colnames(Z_hat), ".z$", "")
  SS <- X %>%
    select(ends_with(".ss")) %>%
    as.matrix()
  nmss <- str_replace(colnames(SS), ".ss$", "")
  o <- match(nmsz, nmss)
  SS <- SS[, o]
  N <- apply(SS, 2, median)
  R <- R_ldsc(Z_hat = Z_hat,
              ldscores = X$L2,
              ld_size = M,
              N = N,
              return_gencov = TRUE,
              return_cor = TRUE,
              make_well_conditioned = TRUE)
  colnames(R$Re) <- rownames(R$Re) <- nmsz
  saveRDS(R$Re, file=paste0(out_dir,prefix,".R_est_ldsc.RDS")) # Double check this
}

#' @title Estimate genetic correlation matrix by LDSC method
#' @param beta_files Paths of combined GWAS data without LD pruning
#' @param ld_files Paths of reference LD score files
#' @param m_files Paths of reference M files
#' @returns A list contain sample overlap matrix (Re) and genetic correlation matrix (Rg)
#'
#' @import ieugwasr
#' @import dplyr
#' @import purrr
#' @import readr
#' @import bigsnpr
#' @import stringr
#' @import sumstatFactors
#' @export
ldsc_full<-function(beta_files, ld_files, m_files){
  ld <- purrr::map_dfr(1:22, function(c){
    read_table(ld_files[c])
  })
  M <- purrr:::map(1:22, function(c){
    read_lines(m_files[c])
  }) %>% unlist() %>% as.numeric() %>% sum()
  X <- map_dfr(beta_files, function(f){
    readRDS(f) %>%
      rename(SNP = snp) %>%
      inner_join(., ld)})
  Z_hat <- X %>%
    select(ends_with(".z")) %>%
    as.matrix()
  nms <- str_replace(colnames(Z_hat), ".z$", "")
  SS <- X %>%
    select(ends_with(".ss")) %>%
    as.matrix()
  N <- apply(SS, 2, median, na.rm=TRUE)

  R <- mrScan::R_ldsc(Z_hat = Z_hat,
              ldscores = X$L2,
              ld_size = M,
              N = N,
              return_gencov = TRUE,
              return_cor = TRUE,
              make_well_conditioned = TRUE)
  Re <- abs(R$Re)
  Rg <- abs(R$Rg)
  colnames(Re) <- rownames(Re) <- colnames(Rg) <- rownames(Rg) <- nms
  return(list(Re = Re, Rg = Rg))
}

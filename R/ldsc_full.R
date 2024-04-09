#' @title Estimate genetic correlation matrix by LDSC method
#' @param dat Combined GWAS data without LD pruning. You have two options:
#' 1. One dataframe of combined data
#' 2. A vector for paths of combined data in each chromosome. It should have 22 elements.
#' eg: c("beta.1.RDS","beta.2.RDS",....,"beta.22.RDS")
#' @param ld_files Paths of reference LD score files
#' @param m_files Paths of reference M files
#' @returns A list contain sample overlap matrix (Re) and genetic correlation matrix (Rg)
#'
#' @import ieugwasr
#' @import dplyr
#' @import readr
#' @import bigsnpr
#' @import stringr
#' @import sumstatFactors
#' @importFrom purrr map_dfr map
#' @export
ldsc_full<-function(dat, ld_files, m_files){
  ld <- purrr::map_dfr(1:22, function(c){
    read_table(ld_files[c])
  })
  M <- purrr:::map(1:22, function(c){
    read_lines(m_files[c])
  }) %>% unlist() %>% as.numeric() %>% sum()
  if(is.vector(dat) & length(dat) == 22){
    X <- map_dfr(dat, function(f){
      readRDS(f) %>%
        rename(SNP = snp) %>%
        inner_join(., ld)})
  }else if(is.data.frame(dat)){
    X <- dat %>%
      rename(SNP = snp) %>%
      inner_join(., ld)
  }
  Z_hat <- X %>%
    select(ends_with(".z")) %>%
    as.matrix()
  nms <- str_replace(colnames(Z_hat), ".z$", "")
  SS <- X %>%
    select(ends_with(".ss")) %>%
    as.matrix()
  N <- apply(SS, 2, median, na.rm=TRUE)
  if(class(N) == "integer"){
    N_int <- N
    N <- as.numeric(N_int)
    names(N) <- names(N_int)
  }
  R <- R_ldsc(Z_hat = Z_hat,
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

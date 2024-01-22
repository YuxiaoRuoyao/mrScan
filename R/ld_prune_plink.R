#' @title LD clumping for harmonized data
#' @param X Paths of merged GWAS summary data by chromosomes
#' @param r2_thresh Clumping r2 cut off. Default = 0.001
#' @param clump_kb Clumping distance cut off. Default = 10000
#' @param type LD clumping prioritization. Either pvalue or rank. Default = "pvalue"
#' @param pthresh pvalue threshold. Default = 1
#' @param ref_path Path for the LD reference panel.
#' @returns One dataframe per chromosome with columns for SNP info
#'
#' @import ieugwasr
#' @import dplyr
#' @export
ld_prune_plink<-function(X,r2_thresh=0.001,clump_kb=10000,type="pvalue",pthresh=1,
                         ref_path){
  miss <- X %>%
    select(ends_with(".z")) %>%
    is.na(.) %>%
    rowSums(.)
  complete_ix <- which(miss == 0)
  X <- X[complete_ix,]

  pmat <- X %>%
    select(ends_with(".p")) %>%
    as.matrix()
  Z_hat <- X %>%
    select(ends_with(".z")) %>%
    as.matrix()
  if(type == "pvalue"){
    myp <- suppressWarnings(apply(as.matrix(pmat[,-1]), 1, function(x){min(x, na.rm=TRUE)}))
  }else if(type == "rank"){
    Z_rank <- apply(Z_hat,2,function(x){rank(x,na.last = "keep")})
    min_rank <- apply(Z_rank, 1, function(x){min(x, na.rm = T)})
    myp <- min_rank/max(min_rank)
  }
  dat <- data.frame(rsid = X$snp, pval = myp)
  dat_clump <- ld_clump(dat = dat,
                        clump_r2 = r2_thresh,
                        clump_p = pthresh,
                        clump_kb = clump_kb,
                        plink_bin = genetics.binaRies::get_plink_binary(),
                        bfile = ref_path)
  ix <- which(X$snp %in% dat_clump$rsid)
  X <- X[ix,]
  return(X)
}

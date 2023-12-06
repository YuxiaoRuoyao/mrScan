#' @title LD clumping for harmonized data
#' @param r2_thresh Clumping r2 cut off. Default = 0.001
#' @param clump_kb Clumping distance cut off. Default = 10000
#' @param p LD clumping prioritization. Either pvalue or rank. Default = "pvalue"
#' @param pthresh pvalue threshold. Default = 1
#' @param zmat_dir Directory path of merged data. Default is the current work directory.
#' @param ref_path Path for the LD reference panel.
#' @param out_dir Output data path. Default is in the current work directory.
#' @param prefix Name prefix for the output. Default = NULL
#' @returns Save one dataframe per chromosome with columns for SNP info
#'
#' @import ieugwasr
#' @import dplyr
#' @export
ld_prune_chrom_plink<-function(r2_thresh=0.001,clump_kb=10000,p="pvalue",pthresh=1,
                               zmat_dir=NULL,prefix=NULL,ref_path,
                               out_dir=NULL){
  for(c in seq(1:22)){
    X <- readRDS(paste0(zmat_dir,prefix,".zmat.",c,".RDS"))
    if(!p %in% c("pvalue", "rank")){
      stop("Unknown prioritization option.\n")
    }
    Z_hat <- X %>%
      select(ends_with(".z")) %>%
      as.matrix()

    if(p == "pvalue"){
      zmax <- apply(abs(Z_hat), 1, function(x){max(x, na.rm=T)})
      myp <- 2*pnorm(-abs(zmax))
    }else if(p == "rank"){
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
    saveRDS(X, file=paste0(out_dir,prefix,".zmat.ldpruned.",c,".RDS"))
  }
}

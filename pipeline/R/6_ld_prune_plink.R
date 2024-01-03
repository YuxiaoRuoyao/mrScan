library(ieugwasr)
library(dplyr)

X <- readRDS(snakemake@input[["beta"]])
r2_thresh <- as.numeric(snakemake@params[["r2_thresh"]])
clump_kb <- snakemake@params[["clump_kb"]]
ref_path  <- snakemake@params[["ref_path"]]
type <- snakemake@params[["ld_prioritization"]]
pthresh <- as.numeric(snakemake@params[["pthresh"]])
out <- snakemake@output[["out"]]


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
  myp <- suppressWarnings(apply(pmat[,-1], 1, function(x){min(x, na.rm=TRUE)}))
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

saveRDS(X, file=out)

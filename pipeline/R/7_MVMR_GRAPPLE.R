library(dplyr)
library(GRAPPLE)
library(purrr)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
R <- readRDS(snakemake@input[["R"]])
R_type <- snakemake@params[["R_type"]]
out <- snakemake@output[["out"]]

if(R_type == "pval"){
  R_matrix <- as.matrix(R)
}else if(R_type == "ldsc"){
  R_matrix <- as.matrix(R$Re)
}

X <- purrr::map_dfr(beta_files, readRDS)
beta_hat <- X %>% select(ends_with(".beta"))
se <- X %>% select(ends_with(".se"))
p <- X %>% select(ends_with(".p"))
nms <- stringr::str_replace(names(beta_hat), ".beta", "")
names(beta_hat)<-names(se)<-names(p)<-nms
o <- match(colnames(R_matrix), nms)
beta_hat <- data.frame(beta_hat[, o],check.names = F)
se <- data.frame(se[, o],check.names = F)
i <- ncol(beta_hat)
grapple_dat <- data.frame(cbind(beta_hat, se))
names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(i-1)),
                        "se_out", paste0("se_exp", 1:(i-1)))
grapple_dat$selection_pvals <- apply(p[,-1, drop = F],1, min)
ix <- which(grapple_dat$selection_pvals < pval_threshold)
if(length(ix) > 5000){
  test_dat<-grapple_dat
  test_dat$index<-seq(1:nrow(grapple_dat))
  select_dat<-head(test_dat[order(test_dat$selection_pvals),],5000)
  ixn<-select_dat$index
  grapple_dat <- grapple_dat[ixn,]
}
res <- grappleRobustEst(data = grapple_dat,
                        plot.it =FALSE,
                        p.thres = pval_threshold,
                        cor.mat = R_matrix,
                        niter = 100000)
res.summary <- data.frame(exposure=colnames(beta_hat)[-1],
                          b=res$beta.hat,
                          se=sqrt(diag(res$beta.var)),
                          pvalue=res$beta.p.value,
                          method = "GRAPPLE")
saveRDS(res.summary,file = out)


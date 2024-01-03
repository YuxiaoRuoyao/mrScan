library(dplyr)
library(MRBEE)
library(purrr)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
pleio_p_thresh <- as.numeric(snakemake@params[["pleio_p_thresh"]])
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
pmin <- apply(p[,-1, drop = F], 1, min)
ix <- which(pmin < pval_threshold)
bT <- list(R = R_matrix, Ncor = Inf,
           EstHarm = beta_hat[ix,],
           SEHarm =  se[ix,])
pD <- prepData(bT,verbose =FALSE)
fit <- MRBEE.IMRP(pD, PleioPThreshold = pleio_p_thresh)
res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                          b = fit$CausalEstimates[-1],
                          se = sqrt(diag(fit$VCovCausalEstimates))[-1])
res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
res.summary$method <- "MRBEE"

saveRDS(res.summary,file = out)

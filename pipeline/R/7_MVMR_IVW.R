library(dplyr)
library(TwoSampleMR)
library(purrr)

beta_files <- unlist(snakemake@input[["beta"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
out <- snakemake@output[["out"]]

X <- purrr::map_dfr(beta_files, readRDS)
beta_hat <- X %>% select(ends_with(".beta"))
se <- X %>% select(ends_with(".se"))
p <- X %>% select(ends_with(".p"))
nms <- stringr::str_replace(names(beta_hat), ".beta", "")
names(beta_hat)<-names(se)<-names(p)<-nms
pmin <- apply(p[,-1, drop = F], 1, min)
ix <- which(pmin < pval_threshold)
i <- ncol(beta_hat)
hdat <-  list(exposure_beta = as.matrix(beta_hat[ix, 2:i]),
              exposure_pval = as.matrix(p[ix, 2:i]),
              exposure_se = as.matrix(se[ix,2:i]),
              outcome_beta = data.frame(beta_hat)[ix,1],
              outcome_pval = data.frame(p)[ix,1],
              outcome_se = data.frame(se)[ix,1],
              expname = data.frame(id.exposure = nms, exposure = nms),
              outname = data.frame(id.outcome = nms[1], outcome = nms[1]))
res_F <- mv_multiple(hdat)$result
res_T <- mv_multiple(hdat,instrument_specific = TRUE)$result
res_F$method <- "IVW"
res_T$method <- "IVW_instrument_specific"

saveRDS(list(IVW_F = res_F, IVW_T = res_T), file = out)

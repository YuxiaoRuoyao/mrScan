library(dplyr)
library(purrr)
library(ebnm)
library(esmr)

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
beta_hat <- X %>% dplyr::select(ends_with(".beta"))
se <- X %>% dplyr::select(ends_with(".se"))
p <- X %>% dplyr::select(ends_with(".p"))
nms <- stringr::str_replace(names(beta_hat), ".beta", "")
names(beta_hat)<-names(se)<-names(p)<-nms
o <- match(colnames(R_matrix), nms)
beta_hat <- data.frame(beta_hat[, o],check.names = F)
se <- data.frame(se[, o],check.names = F)
i <- ncol(beta_hat)
fit <- esmr(beta_hat_Y <- beta_hat[,1],
            se_Y <- se[,1],
            beta_hat_X <- beta_hat[,2:i],
            se_X <- se[, 2:i],
            R = R_matrix,
            augment_G = TRUE,
            g_type = "gfa",
            ix1 = "pval-5e-8",
            ix0 = FALSE,
            lfsr_thresh = 1)

res.summary <- data.frame(exposure = colnames(beta_hat)[-1],
                          b = fit$beta$beta_m,
                          se = fit$beta$beta_s)
res.summary$pvalue <- with(res.summary, 2*pnorm(-abs(b/se)))
res.summary$method <- paste0("ESMR_",pval_threshold)

saveRDS(res.summary,file = out)

library(dplyr)
library(purrr)
library(stringr)
library(sumstatFactors)

beta_files <- unlist(snakemake@input[["beta"]])
p_thresh <- as.numeric(snakemake@params[["p_thresh"]])
out <- snakemake@output[["out"]]

X <- purrr::map_dfr(beta_files, readRDS)
beta_hat <- X %>%
  select(ends_with(".beta")) %>%
  as.matrix()
se_hat <- X %>%
  select(ends_with(".se")) %>%
  as.matrix()
nms <- colnames(beta_hat) %>% stringr::str_replace(".beta$", "")
Rpt <- R_pt(B_hat = beta_hat,
            S_hat = se_hat,
            p_val_thresh = p_thresh,
            return_cor = TRUE,
            make_well_conditioned = TRUE
)
colnames(Rpt) <- rownames(Rpt) <- nms
saveRDS(Rpt, file=out)

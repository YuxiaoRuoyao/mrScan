library(hdme)
library(dplyr)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
radius_type <- snakemake@params[["radius_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
maxits <- as.numeric(snakemake@params[["maxits"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
mvdat_y <- mvdat$mvdat_y

if(radius_type == "min"){
  radius <- 'radius_min'
}else{
  radius <- 'radius_1se'
}
ids <- colnames(mvdat_y$exposure_beta)
ids <- ids[!ids %in% id_exposure]
ss <- df_info %>% filter(id %in% ids) %>%
  arrange(match(id, ids)) %>% pull(sample_size)
set.seed(seed)
cv_corrected_y <- cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                     y = mvdat_y$outcome_beta,
                                     sigmaUU = diag(1/ss))
corrected_y <- corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                     y = mvdat_y$outcome_beta,
                                     sigmaUU = diag(1/ss),
                                     radii = cv_corrected_y[[radius]],
                                     maxits = maxits)
coef_y <-coef(corrected_y)
id.select <- ids[coef_y$coefficient]
df_info[df_info$id %in% id.select,"status"] <- "Select by Corrected Lasso"

saveRDS(list(id.list=id.select,trait.info=df_info),file = out)

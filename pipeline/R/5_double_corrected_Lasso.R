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
mvdat_x <- mvdat$mvdat_x
mvdat_y <- mvdat$mvdat_y

if(radius_type == "min"){
  radius <- 'radius_min'
}else{
  radius <- 'radius_1se'
}
ids_x <- colnames(mvdat_x$exposure_beta)
ss_x <- df_info %>% filter(id %in% ids_x) %>%
  arrange(match(id, ids_x)) %>% pull(sample_size)
set.seed(seed)
cv_corrected_x <- cv_corrected_lasso(W = mvdat_x$exposure_beta,
                                     y = mvdat_x$outcome_beta,
                                     sigmaUU = diag(1/ss_x))
corrected_x <- corrected_lasso(W = mvdat_x$exposure_beta,
                               y = mvdat_x$outcome_beta,
                               sigmaUU = diag(1/ss_x),
                               radii = cv_corrected_x[[radius]],
                               maxits = maxits)
coef_x <-coef(corrected_x)
id_x <- ids_x[coef_x$coefficient]

ids_y <- colnames(mvdat_y$exposure_beta)
ids_y <- ids_y[!ids_y %in% id_exposure]
ss_y <- df_info %>% filter(id %in% ids_y) %>%
  arrange(match(id, ids_y)) %>% pull(sample_size)
set.seed(seed)
cv_corrected_y <- cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                     y = mvdat_y$outcome_beta,
                                     sigmaUU = diag(1/ss_y))
corrected_y <- corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                               y = mvdat_y$outcome_beta,
                               sigmaUU = diag(1/ss_y),
                               radii = cv_corrected_y[[radius]],
                               maxits = maxits)
coef_y <-coef(corrected_y)
id_y <- ids_y[coef_y$coefficient]
id.select <- union(id_x,id_y)
df_info[df_info$id %in% id.select,"status"] <- "Select by Double Corrected Lasso"

saveRDS(list(id.list=id.select,trait.info=df_info),file = out)

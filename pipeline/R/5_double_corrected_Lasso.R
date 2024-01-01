library(TwoSampleMR)
library(hdme)
library(dplyr)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
find_proxies <- as.logical(snakemake@params[["find_proxies"]])
pop <- snakemake@params[["population"]]
harmonise_strictness <- as.numeric(snakemake@params[["harmonise_strictness"]])
radius_type <- snakemake@params[["radius_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
maxits <- as.numeric(snakemake@params[["maxits"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info

if(radius_type == "min"){
  radius <- 'radius_min'
}else{
  radius <- 'radius_1se'
}
# X ~ Z
inst_x <- TwoSampleMR::mv_extract_exposures(id.list,clump_r2 = r2, clump_kb = kb,
                                            harmonise_strictness = harmonise_strictness,find_proxies = find_proxies,
                                            pval_threshold = pval_threshold, pop = pop)
out_x <- TwoSampleMR::extract_outcome_data(inst_x$SNP, id_exposure)
mvdat_x <- mv_harmonise_data(inst_x, out_x)
ids_x <- colnames(mvdat_x$exposure_beta)
ss_x <- df_info %>% filter(id %in% ids_x) %>%
  arrange(match(id, ids_x)) %>% pull(sample_size)
set.seed(seed)
cv_corrected_x <- hdme::cv_corrected_lasso(W = mvdat_x$exposure_beta,
                                           y = mvdat_x$outcome_beta,
                                           sigmaUU = diag(1/ss_x))
corrected_x <- hdme::corrected_lasso(W = mvdat_x$exposure_beta,
                                     y = mvdat_x$outcome_beta,
                                     sigmaUU = diag(1/ss_x),
                                     radii = cv_corrected_x[[radius]],
                                     maxits = maxits)
coef_x <-coef(corrected_x)
id_x <- ids_x[coef_x$coefficient]
# Y ~ Z
inst_y <- TwoSampleMR::mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                            harmonise_strictness = harmonise_strictness,
                                            find_proxies = find_proxies,
                                            pval_threshold = pval_threshold, pop = pop)
out_y <- TwoSampleMR::extract_outcome_data(inst_y$SNP, id_outcome)
mvdat_y <- mv_harmonise_data(inst_y, out_y)
ids_y <- colnames(mvdat_y$exposure_beta)
ids_y <- ids_y[!ids_y %in% id_exposure]
ss_y <- df_info %>% filter(id %in% ids_y) %>%
  arrange(match(id, ids_y)) %>% pull(sample_size)
set.seed(seed)
cv_corrected_y <- hdme::cv_corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                           y = mvdat_y$outcome_beta,
                                           sigmaUU = diag(1/ss_y))
corrected_y <- hdme::corrected_lasso(W = mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                                     y = mvdat_y$outcome_beta,
                                     sigmaUU = diag(1/ss_y),
                                     radii = cv_corrected_y[[radius]],
                                     maxits = maxits)
coef_y <-coef(corrected_y)
id_y <- ids_y[coef_y$coefficient]
id.select <- union(id_x,id_y)
df_info[df_info$id %in% id.select,"status"] <- "Select by Double Corrected Lasso"

saveRDS(list(id.list=id.select,trait.info=df_info),file = out)

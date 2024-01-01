library(TwoSampleMR)
library(glmnet)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
find_proxies <- as.logical(snakemake@params[["find_proxies"]])
pop <- snakemake@params[["population"]]
harmonise_strictness <- as.numeric(snakemake@params[["harmonise_strictness"]])
lambda_type <- snakemake@params[["lambda_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info

if(lambda_type == "min"){
  penalty <- 'lambda.min'
}else{
  penalty <- 'lambda.1se'
}
# X ~ Z
inst_x <- mv_extract_exposures(id.list,clump_r2 = r2, clump_kb = kb,
                                            harmonise_strictness = harmonise_strictness,find_proxies = find_proxies,
                                            pval_threshold = pval_threshold, pop = pop)
out_x <- extract_outcome_data(inst_x$SNP, id_exposure)
mvdat_x <- mv_harmonise_data(inst_x, out_x)
set.seed(seed)
cv_model_x <- cv.glmnet(x=mvdat_x$exposure_beta, y=mvdat_x$outcome_beta,
                                alpha = 1)
best_model_x <- glmnet(x=mvdat_x$exposure_beta,
                               y=mvdat_x$outcome_beta, alpha = 1,
                               lambda = cv_model_x[[penalty]])
A_x<-rownames(coef(best_model_x))[which(coef(best_model_x)[,1]!=0)][-1]
# Y ~ Z
inst_y <- mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                            harmonise_strictness = harmonise_strictness,
                                            find_proxies = find_proxies,
                                            pval_threshold = pval_threshold, pop = pop)
out_y <- extract_outcome_data(inst_y$SNP, id_outcome)
mvdat_y <- mv_harmonise_data(inst_y, out_y)
set.seed(seed)
cv_model_y <- cv.glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                        y=mvdat_y$outcome_beta, alpha = 1)
best_model_y <- glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                       y=mvdat_y$outcome_beta, alpha = 1,
                       lambda = cv_model_y[[penalty]])
A_y<-rownames(coef(best_model_y))[which(coef(best_model_y)[,1]!=0)][-1]
A<-union(A_x,A_y)
df_info[df_info$id %in% A,"status"] <- "Select by Double Lasso"

saveRDS(list(id.list=A,trait.info=df_info),file = out)

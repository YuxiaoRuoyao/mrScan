library(TwoSampleMR)
library(glmnet)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
mvdat <- readRDS(snakemake@input[["instruments"]])
lambda_type <- snakemake@params[["lambda_type"]]
seed <- as.numeric(snakemake@params[["seed"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info
mvdat_y <- mvdat$mvdat_y

if(lambda_type == "min"){
  penalty <- 'lambda.min'
}else{
  penalty <- 'lambda.1se'
}
set.seed(seed)
# not penalize X
cv_model_y <- cv.glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                        y=mvdat_y$outcome_beta, alpha = 1)
best_model_y <- glmnet(x=mvdat_y$exposure_beta[,-which(colnames(mvdat_y$exposure_beta)==id_exposure)],
                       y=mvdat_y$outcome_beta, alpha = 1,
                       lambda = cv_model_y[[penalty]])
A_y<-rownames(coef(best_model_y))[which(coef(best_model_y)[,1]!=0)][-1]
df_info[df_info$id %in% A_y,"status"] <- "Select by Classic Lasso"
saveRDS(list(id.list=A_y,trait.info=df_info),file = out)

library(TwoSampleMR)
library(stringr)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
find_proxies <- as.logical(snakemake@params[["find_proxies"]])
pop <- snakemake@params[["population"]]
harmonise_strictness <- as.numeric(snakemake@params[["harmonise_strictness"]])
method <- snakemake@params[["method"]]
out <- snakemake@output[["out"]]

id.list <- res$id.list
df_info <- res$trait.info

inst_y <- mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                                            harmonise_strictness = harmonise_strictness,
                                            find_proxies = find_proxies,
                                            pval_threshold = pval_threshold, pop = pop)
out_y <- extract_outcome_data(inst_y$SNP, id_outcome)
mvdat_y <- mv_harmonise_data(inst_y, out_y)
dat_all <- data.frame(mvdat_y$exposure_beta)
names(dat_all) <- colnames(mvdat_y$exposure_beta)
dat_all$y <- mvdat_y$outcome_beta
w <- 1/mvdat_y$outcome_se^2
if(method == "forward"){
  fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
  f1 <- stats::lm(fo1, data = dat_all, weights = w)
  fall <- stats::lm(y ~ ., data = dat_all, weights = w)
  fstep <- stats::step(f1, scope = formula(fall), direction = "forward",trace = FALSE)
  ii <- summary(fstep)$coefficients %>% rownames()
  ii <- ii[-1]
  id.select <- str_replace_all(ii, stringr::fixed("`"), "")
}else if(method == "backward"){
  fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
  f1 <- stats::lm(fo1, data = dat_all, weights = w)
  fall <- stats::lm(y ~ ., data = dat_all, weights = w)
  fback <- stats::step(fall,scope=list(upper=formula(fall),lower=formula(f1)),
                       direction="backward",trace=FALSE)
  ii <- summary(fback)$coefficients %>% rownames()
  ii <- ii[-1]
  id.select <- str_replace_all(ii, stringr::fixed("`"), "")
}else if(method == "both"){
  fo1 <- as.formula(paste0("y ~ `", id_exposure,"`"))
  f1 <- stats::lm(fo1, data = dat_all, weights = w)
  fall <- stats::lm(y ~ ., data = dat_all, weights = w)
  fboth <- stats::step(f1,scope=list(upper=formula(fall),lower=formula(f1)),
                       direction="both",trace=FALSE)
  ii <- summary(fboth)$coefficients %>% rownames()
  ii <- ii[-1]
  id.select <- str_replace_all(ii, stringr::fixed("`"), "")
}
id.select <- id.select[!id.select %in% id_exposure]
df_info[df_info$id %in% id.select,"status"] <- paste0("Select by stepwise ",method)

saveRDS(list(id.list=id.select,trait.info=df_info),file = out)

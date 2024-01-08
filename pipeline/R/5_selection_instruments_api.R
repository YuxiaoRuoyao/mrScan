library(TwoSampleMR)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
res <- readRDS(snakemake@input[["file"]])
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
find_proxies <- as.logical(snakemake@params[["find_proxies"]])
pop <- snakemake@params[["population"]]
harmonise_strictness <- as.numeric(snakemake@params[["harmonise_strictness"]])
out <- snakemake@output[["out"]]

id.list <- res$id.list

# X ~ Z
inst_x <- mv_extract_exposures(id.list,clump_r2 = r2, clump_kb = kb,
                               harmonise_strictness = harmonise_strictness,find_proxies = find_proxies,
                               pval_threshold = pval_threshold, pop = pop)
out_x <- extract_outcome_data(inst_x$SNP, id_exposure)
mvdat_x <- mv_harmonise_data(inst_x, out_x)
# Y ~ Z
inst_y <- mv_extract_exposures(c(id_exposure,id.list),clump_r2 = r2, clump_kb = kb,
                               harmonise_strictness = harmonise_strictness,
                               find_proxies = find_proxies,
                               pval_threshold = pval_threshold, pop = pop)
out_y <- extract_outcome_data(inst_y$SNP, id_outcome)
mvdat_y <- mv_harmonise_data(inst_y, out_y)

saveRDS(list(mvdat_x=mvdat_x,mvdat_y=mvdat_y),file = out)

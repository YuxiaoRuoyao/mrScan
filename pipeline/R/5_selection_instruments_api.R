library(TwoSampleMR)
library(mrScan)

id_exposure <- snakemake@params[["id_exposure"]]
id_outcome <- snakemake@params[["id_outcome"]]
id.list <- read.csv(snakemake@input[["id_list"]])$id
r2 <- as.numeric(snakemake@params[["r2_thresh"]])
kb <- as.numeric(snakemake@params[["clump_kb"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
find_proxies <- as.logical(snakemake@params[["find_proxies"]])
pop <- snakemake@params[["population"]]
harmonise_strictness <- as.numeric(snakemake@params[["harmonise_strictness"]])
out <- snakemake@output[["out"]]

res_inst_api <- select_instruments_api(id.list = id.list,id_exposure = id_exposure,
                                       id_outcome = id_outcome,r2 = r2, kb = kb,
                                       pval_threshold = pval_threshold,
                                       find_proxies = find_proxies,pop = pop,
                                       harmonise_strictness = harmonise_strictness)
saveRDS(res_inst_api,file = out)

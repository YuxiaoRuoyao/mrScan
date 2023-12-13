library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(purrr)
inst_files <- unlist(snakemake@input[["inst_files"]])
trait <- snakemake@params[["trait"]]
out <- snakemake@output[["out"]]

all_inst <- map(inst_files, function(f){
  readRDS(f)
}) %>% unlist() %>% unique()
chunk_inst <- split(all_inst, ceiling(seq_along(all_inst)/20))
f <- map_dfr(chunk_inst, function(x) ieugwasr::associations(x,trait) %>%
                      distinct(rsid, .keep_all = TRUE))
new_f <- data.frame(SNP = all_inst) %>%
  left_join(f[,c("rsid","beta")],by=c("SNP" = "rsid"))
colnames(new_f)[2] <- trait
saveRDS(new_f,file=out)

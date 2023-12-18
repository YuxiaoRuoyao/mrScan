library(ieugwasr)
library(dplyr)
source("R/helpers.R")
res_initial <- snakemake@input[["file"]]
id_exposure <- snakemake@params[["id_exposure"]]
nsnp_cutoff <- as.numeric(snakemake@params[["nsnp_cutoff"]])
pop <- snakemake@params[["population"]]
sex <- snakemake@params[["sex"]]
out_id_list <- snakemake@output[["id_list"]]
out_trait_info <- snakemake@output[["trait_info"]]

dat <- readRDS(res_initial)$trait.info

na.SNP.trait<-filter(dat, is.na(nsnp))$id
na.sex.trait<-filter(dat,sex=="NA")$id
new.dat<-filter(dat, nsnp > nsnp_cutoff & sex == sex & population == pop)
id.list <- new.dat$id
id.list<-unique(c(id.list,na.SNP.trait,na.sex.trait))
new.dat<-dat %>% filter(id %in% id.list) %>% filter(population == pop)
id.list<-new.dat$id
dat <- dat %>%
  mutate(status = if_else(id %in% id.list, "select after QC", "delete in QC"))
# delete menarche trait with sex males and females
dat <- dat %>%
  mutate(status = ifelse(grepl('menarche', trait) == TRUE & sex == "Males and Females",
                                "delete in QC",status))
# delete trait-adjusted traits
dat <- dat %>%
  mutate(status = ifelse(grepl('adjust', trait) == TRUE,
                                "delete in QC",status))
id.qc <- dat %>% filter(status == "select after QC") %>% pull(id)

write.csv(data.frame(id = id.qc),file = out_id_list,row.names = F)
write.csv(dat,file = out_trait_info,row.names = F)

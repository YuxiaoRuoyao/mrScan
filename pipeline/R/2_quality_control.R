library(ieugwasr)
library(dplyr)
source("R/helpers.R")
res_initial <- snakemake@input
id_exposure <- snakemake@params[["id_exposure"]]
nsnp_cutoff <- as.numeric(snakemake@params[["id_exposure"]])
pop <- snakemake@params[["population"]]
sex <- snakemake@params[["sex"]]
R2_cutoff <- as.numeric(snakemake@params[["R2_cutoff"]])
out_id_list <- snakemake@output[["id_list"]]
out_trait_info <- snakemake@output[["trait_info"]]

dat <- readRDS(res_initial)$trait.info

na.SNP.trait<-filter(dat, is.na(nsnp))$id
na.sex.trait<-filter(dat,sex=="NA")$id
new.dat<-filter(new.dat, nsnp > nsnp_cutoff & sex == sex & population == pop)
id.list <- new.dat$id
id.list<-unique(c(id.list,na.SNP.trait,na.sex.trait))
new.dat<-info %>% filter(id %in% id.list) %>% filter(population == pop)
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
# delete high correlation traits with X
res_cor_X <- calculate_cor(ids1 = id_exposure, ids2 = id.list)
trait_cor_X <- filter(res_cor_X, abs(cor) > R2_cutoff) %>%
  with(., c(id1, id2) ) %>% unique()
dat <- dat %>%
  mutate(status = ifelse(id %in% trait_cor_X,
                                "delete since high cor with X",status))
id.qc <- dat %>% filter(status == "select after QC") %>% pull(id)

write.csv(data.frame(id = id.qc),file = out_id_list,row.names = F)
write.csv(dat,file = out_trait_info,row.names = F)

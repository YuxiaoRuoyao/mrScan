#' @title Quality control based on metadata
#' @param id_exposure GWAS ID for main exposure.
#' @param dat metadata for traits, a dataframe contains information of population, nsnp and etc.
#' If users input metadata by themselves, use make_metadata() to convert data format
#' @param nsnp_cutoff threshold of minimum number of SNPs. Default = 1e6
#' @param pop Limit population of traits. Default = "European"
#' @param sex Gender limitation. Default = "Males and Females"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import mrScan
#' @import ieugwasr
#' @import dplyr
#' @export
quality_control <- function(id_exposure,dat,nsnp_cutoff=1e6,pop="European",sex="Males and Females"){
  na.SNP.trait<-filter(dat, is.na(nsnp))$id
  na.sex.trait<-filter(dat,sex=="NA")$id
  new.dat<-filter(dat, population == pop)
  new.dat<-filter(new.dat, nsnp > nsnp_cutoff & sex == sex)
  id.list <- new.dat$id
  id.list<-unique(c(id.list,na.SNP.trait,na.sex.trait))
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
  res_cor_X <- mr_pairs_raps(ids1 = id_exposure, ids2 = id.list)
  trait_cor_X <- filter(res_cor_X$cor, abs(cor) > 0.9) %>%
    with(., c(id1, id2) ) %>% unique()
  dat <- dat %>%
    mutate(status = ifelse(id %in% trait_cor_X,
                           "delete since high cor with X",status))
  id.qc <- dat %>% filter(status == "select after QC") %>% pull(id)
  return(list(id.list=id.qc,trait.info=dat))
}

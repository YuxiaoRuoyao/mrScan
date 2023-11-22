#' @title Quality control based on metadata
#' @param id_exposure GWAS ID for main exposure.
#' @param dat metadata for traits, a dataframe contains information of population, nsnp and etc.
#' If users input metadata by themselves, use make_metadata() to convert data format
#' @param nsnp_cutoff threshold of minimum number of SNPs. Default = 1e6
#' @param pop Limit population of traits. Default = "European"
#' @param sex Gender limitation. Default = "Males and Females"
#' @param R2_cutoff cutoff for duplicated traits of the exposure. Default = 0.9
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom dplyr filter mutate
#' @export
quality_control <- function(id_exposure,dat,nsnp_cutoff=1e6,pop="European",sex="Males and Females",
                            R2_cutoff = 0.9){
  na.SNP.trait<-dplyr::filter(dat, is.na(nsnp))$id
  na.sex.trait<-dplyr::filter(dat,sex=="NA")$id
  new.dat<-dplyr::filter(dat, population == pop)
  new.dat<-dplyr::filter(new.dat, nsnp > nsnp_cutoff & sex == sex)
  id.list <- new.dat$id
  id.list<-unique(c(id.list,na.SNP.trait,na.sex.trait))
  dat <- dat %>%
    mutate(status = if_else(id %in% id.list, "select after QC", "delete in QC"))
  # delete menarche trait with sex males and females
  dat <- dat %>%
    dplyr::mutate(status = ifelse(grepl('menarche', trait) == TRUE & sex == "Males and Females",
                           "delete in QC",status))
  # delete trait-adjusted traits
  dat <- dat %>%
    dplyr::mutate(status = ifelse(grepl('adjust', trait) == TRUE,
                           "delete in QC",status))
  # delete high correlation traits with X
  res_cor_X <- calculate_cor(ids1 = id_exposure, ids2 = id.list)
  trait_cor_X <- dplyr::filter(res_cor_X, abs(cor) > R2_cutoff) %>%
    with(., c(id1, id2) ) %>% unique()
  dat <- dat %>%
    dplyr::mutate(status = ifelse(id %in% trait_cor_X,
                           "delete since high cor with X",status))
  id.qc <- dat %>% dplyr::filter(status == "select after QC") %>% pull(id)
  return(list(id.list=id.qc,trait.info=dat))
}
# add option to let users filter by themselves

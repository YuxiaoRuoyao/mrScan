#' @title Quality control based on metadata
#' @param dat metadata for traits, a dataframe contains information of population, nsnp and etc.
#' If users input metadata by themselves, use make_metadata() to convert data format
#' @param nsnp_cutoff threshold of minimum number of SNPs. Default = 1e6
#' @param pop Limit population of traits. Default = "European"
#' @param gender Gender limitation. Default = "Males and Females"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @export
quality_control <- function(dat,nsnp_cutoff=1e6,pop="European",gender="Males and Females"){
  id.list1 <- dat %>% filter(nsnp > nsnp_cutoff & sex == gender & population == pop) %>% pull(id)
  id.list2 <- dat %>% filter(is.na(nsnp) & sex == gender & population == pop) %>% pull(id)
  id.list3 <- dat %>% filter(is.na(sex) & nsnp > nsnp_cutoff & population == pop) %>% pull(id)
  id.list <- unique(c(id.list1,id.list2,id.list3))
  dat <- dat %>%
    dplyr::mutate(status = if_else(id %in% id.list, "select after QC", "delete in QC"))
  # delete menarche trait with sex males and females
  dat <- dat %>%
    dplyr::mutate(status = ifelse(grepl('menarche', trait) == TRUE & sex == "Males and Females",
                           "delete in QC",status))
  # delete trait-adjusted traits
  dat <- dat %>%
    dplyr::mutate(status = ifelse(grepl('adjust', trait) == TRUE,
                           "delete in QC",status))
  id.qc <- dat %>% dplyr::filter(status == "select after QC") %>% pull(id)
  return(list(id.list=id.qc,trait.info=dat))
}

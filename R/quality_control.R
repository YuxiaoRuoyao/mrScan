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
  id.list3 <- dat %>% filter(is.na(sex) | sex == "NA" & nsnp > nsnp_cutoff & population == pop) %>% pull(id)
  id.list <- unique(c(id.list1,id.list2,id.list3))
  id_exposure <- dat %>% filter(status == "Delete due to it's the main exposure") %>% pull(id)
  id.list <- id.list[!id.list %in% id_exposure]
  dat <- dat %>% mutate(status = case_when(id %in% id.list ~ "select after QC",
                                           id %in% id_exposure  ~ "Delete due to it's the main exposure",
                                           TRUE ~ "delete in QC"))
  # delete menarche trait with sex males and females
  dat <- dat %>%
    mutate(status = ifelse(grepl('menarche', trait, ignore.case = TRUE), "delete in QC", status))
  # delete trait-adjusted traits
  dat <- dat %>%
    mutate(status = ifelse(grepl('adjust', trait) == TRUE,
                           "delete in QC",status))
  id.qc <- dat %>% filter(status == "select after QC") %>% pull(id)
  return(list(id.list=id.qc,trait.info=dat))
}

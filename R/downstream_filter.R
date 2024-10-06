#' @title Filter downstream traits
#' @description Use Bidirectional MR or MVMR results to filter out downstream traits.
#' It's based on FDR-adjusted pvalue on four directions. Select either upstream traits
#' for exposure or outcome direction (p_ZtoX_adj < p | p_ZtoY_adj < p) and
#' exclude all downstream traits for both directions (p_XtoZ_adj < p | p_YtoZ_adj < p).
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param res Bidirection MR results for selected traits
#' @param p p-value cutoff for adjusted pval to select upstream or downstream traits. Default=0.05
#' @param MR_method Bidirection MR or MVMR methods. Options: MR_IVW, MR_GRAPPLE,
#' MR_MRBEE, MVMR_IVW, MVMR_GRAPPLE, MVMR_MRBEE
#' @param extra_traits trait ID which is adjusted for in bidirection MR. Default = "None"
#' @returns A GWAS ID vector, a trait info dataframe, a trait dataframe with four direction
#' estimate and t-test results
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom data.table dcast setDT
#' @export
downstream_filter <- function(id_exposure,id.list,df_info,res,p = 0.05,
                              extra_traits = "None",MR_method){
  df_summary<-data.frame(id=id.list) %>% left_join(df_info[,c("id","trait")],by="id")
  df_YtoZ <- left_join(df_summary,res$mr12,by=c("id"="id.outcome")) %>%
    filter(method == MR_method) %>%
    mutate(exposure = if_else(id.exposure == id_exposure, "X", "Y"))
  df_ZtoY <- left_join(df_summary,res$mr21,by=c("id"="id.exposure")) %>%
    filter(method == MR_method) %>%
    mutate(outcome = if_else(id.outcome == id_exposure, "X", "Y"))
  notconverge_trait <- character(0)
  if(MR_method == "MR_GRAPPLE" | MR_method == "MVMR_GRAPPLE"){
    id_notconverge1 <- df_YtoZ[df_YtoZ$converge == FALSE | df_YtoZ$se > 1,"id"]
    id_notconverge2 <- df_ZtoY[df_ZtoY$converge == FALSE | df_ZtoY$se > 1,"id"]
    notconverge_trait <- intersect(id_notconverge1,id_notconverge2)
    df_info[df_info$id %in% notconverge_trait, "status"] <- paste0("delete due to not converge by ", MR_method)
    df_YtoZ <- df_YtoZ %>% filter(converge == TRUE & se <= 1)
    df_ZtoY <- df_ZtoY %>% filter(converge == TRUE & se <= 1)
  }else{
    id_notconverge1 <- df_YtoZ[df_YtoZ$se > 1,"id"]
    id_notconverge2 <- df_ZtoY[df_ZtoY$se > 1,"id"]
    notconverge_trait <- intersect(id_notconverge1,id_notconverge2)
    df_info[df_info$id %in% notconverge_trait, "status"] <- paste0("delete due to large SE by ", MR_method)
    df_YtoZ <- df_YtoZ %>% filter(se <= 1)
    df_ZtoY <- df_ZtoY %>% filter(se <= 1)
  }
  setDT(df_YtoZ)
  data_wide1 <- dcast(df_YtoZ,id ~ exposure, value.var=c("b","se","pvalue")) %>%
    select_if(~sum(!is.na(.)) > 0)
  colnames(data_wide1)<-c("id","b_XtoZ","b_YtoZ","se_XtoZ","se_YtoZ","p_XtoZ","p_YtoZ")
  setDT(df_ZtoY)
  data_wide2 <- dcast(df_ZtoY,id ~ outcome, value.var=c("b","se","pvalue")) %>%
    select_if(~sum(!is.na(.)) > 0)
  colnames(data_wide2)<-c("id","b_ZtoX","b_ZtoY","se_ZtoX","se_ZtoY","p_ZtoX","p_ZtoY")
  df_final<-full_join(data_wide1,data_wide2,by="id") %>%
    mutate(across(c("p_ZtoX","p_XtoZ","p_ZtoY","p_YtoZ"),
                  ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_adj"))
  sig_trait <- df_final %>% filter(p_ZtoX_adj < p | p_ZtoY_adj < p) %>% pull(id)
  down_trait <- df_final %>% filter(p_XtoZ_adj < p | p_YtoZ_adj < p) %>% pull(id)
  select_trait <- sig_trait[!sig_trait %in% c(down_trait,notconverge_trait)]
  df_info[df_info$id %in% select_trait,"status"] <- paste0("select after downstream filtering by ",MR_method)
  filter.trait <- id.list[!id.list %in% c(select_trait,notconverge_trait)]
  df_info[df_info$id %in% filter.trait,"status"] <- paste0("delete in downstream filtering by ",MR_method)
  if(extra_traits != "None"){
    select_trait <- c(select_trait,extra_traits)
    df_info[df_info$id %in% extra_traits,"status"] <- "select after downstream filtering"
  }
  return(list(id.list=select_trait,trait.info=df_info,df_bidirection = df_final))
}

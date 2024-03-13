#' @title Filter downstream traits
#' @param id_exposure GWAS ID of the main exposure
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param res Bidirection MR results for selected traits
#' @param sig_level One-sided t-test significant level. Default=0.05
#' @param MR_method Bidirection MR or MVMR methods. Options: MR_IVW, MR_GRAPPLE,
#' MR_MRBEE, MVMR_IVW, MVMR_GRAPPLE, MVMR_MRBEE
#' @returns A GWAS ID vector, a trait info dataframe, a trait dataframe with four direction
#' estimate and t-test results
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom data.table dcast setDT
#' @export
downstream_filter <- function(id_exposure,id.list,df_info,res,sig_level = 0.05,MR_method){
  df_summary<-data.frame(id=id.list) %>% left_join(df_info[,c("id","trait")],by="id")
  df_YtoZ <- left_join(df_summary,res$mr12,by=c("id"="id.outcome")) %>%
    filter(method == MR_method) %>%
    mutate(exposure = if_else(id.exposure == id_exposure, "X", "Y"))
  setDT(df_YtoZ)
  data_wide1 <- dcast(df_YtoZ,id ~ exposure, value.var=c("b","se","pvalue")) %>%
    select_if(~sum(!is.na(.)) > 0)
  colnames(data_wide1)<-c("id","b_XtoZ","b_YtoZ","se_XtoZ","se_YtoZ","p_XtoZ","p_YtoZ")
  df_ZtoY <- left_join(df_summary,res$mr21,by=c("id"="id.exposure")) %>%
    filter(method == MR_method) %>%
    mutate(outcome = if_else(id.outcome == id_exposure, "X", "Y"))
  setDT(df_ZtoY)
  data_wide2 <- dcast(df_ZtoY,id ~ outcome, value.var=c("b","se","pvalue")) %>%
    select_if(~sum(!is.na(.)) > 0)
  colnames(data_wide2)<-c("id","b_ZtoX","b_ZtoY","se_ZtoX","se_ZtoY","p_ZtoX","p_ZtoY")
  df_final<-left_join(data_wide1,data_wide2,by="id")
  df_final<- df_final %>%
    mutate(t_Z_X = (abs(b_ZtoX)-abs(b_XtoZ))/sqrt(se_ZtoX^2+se_XtoZ^2),
           t_Z_Y = (abs(b_ZtoY)-abs(b_YtoZ))/sqrt(se_ZtoY^2+se_YtoZ^2)) %>%
    mutate(p_Z_X = pnorm(t_Z_X),p_Z_Y = pnorm(t_Z_Y))
  X_downstream<- df_final %>% filter(p_Z_X < sig_level) %>% pull(id)
  X_sig<- df_final %>% filter(p_Z_X > 1 - sig_level) %>% pull(id)
  Y_downstream<- df_final %>% filter(p_Z_Y < sig_level) %>% pull(id)
  Y_sig<- df_final %>% filter(p_Z_Y > 1 - sig_level) %>% pull(id)
  sig_traits <- union(X_sig,Y_sig)
  downstream_traits <- union(X_downstream,Y_downstream)
  select_trait <- sig_traits[!sig_traits %in% downstream_traits]
  df_info[df_info$id %in% select_trait,"status"] <- paste0("select after downstream filtering by ",MR_method)
  filter.trait <- id.list[!id.list %in% select_trait]
  df_info[df_info$id %in% filter.trait,"status"] <- paste0("delete in downstream filtering by ",MR_method)
  return(list(id.list=select_trait,trait.info=df_info,df_bidirection = df_final))
}

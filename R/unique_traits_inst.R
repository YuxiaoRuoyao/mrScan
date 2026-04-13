#' @title Select unique traits based on pairwise genetic correlation
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param R_matrix Pairwise genetic correlation matrix. Colnames are trait names
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=0.8
#' @param extra_traits trait ID which is adjusted for in bidirection MR. Default = "None"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom dplyr left_join pull arrange
#' @export
unique_traits_inst <- function(id.list,df_info,R_matrix,R2_cutoff=0.8,extra_traits = "None"){
  res_group <- greedy_group(id.list = id.list,R = R_matrix,
                            df_info = df_info,R2_cutoff = R2_cutoff)
  clusters <- res_group$cluster_info
  df_info <- left_join(df_info,clusters,by = c("id" = "id"))
  df_cluster <- left_join(clusters,df_info[,c("id","n_inst")],
                          by=c("id" = "id")) %>% arrange(cluster)
  ids.final <- res_group$selected_traits
  if(extra_traits != "None"){
    ids.final <- df_cluster %>% filter(id %in% extra_traits) %>% 
      bind_rows(df_cluster %>% filter(id %in% ids.final)) %>%
      distinct(cluster, .keep_all = T) %>% arrange(cluster) %>% pull(id)
  }
  filter.trait <- df_cluster$id[!df_cluster$id %in% ids.final]
  df_info[df_info$id %in% filter.trait,"status"] <- "delete due to duplicates"
  df_info[df_info$id %in% ids.final,"status"] <- "select after filtering duplicates"
  return(list(id.list=ids.final,trait.info=df_info))
}

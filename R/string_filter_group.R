#' @title Use string similarity to filter traits
#' @description Keep traits with the highest number of instruments in one cluster. The trait correlation matrix is based on string similarity score and then use greedy cluster.
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param df_inst_counts Dataframe of the number of instruments for each trait.
#' The required columns include `id` for trait ID, `n_inst` for the number of instruments.
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @param extra_traits trait ID which you must want to include and select. Default = "None"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import dplyr
#' @export
string_filter_group <- function(id.list,df_info,df_inst_counts,R2_cutoff = 0.9,
                                extra_traits = "None"){
  res_string <- string_sim(id.list = id.list,df_info = df_info)
  res_group <- greedy_group(id.list = id.list, R = res_string$R_matrix,
                            df_info = df_inst_counts,
                            R2_cutoff = R2_cutoff)
  res_cluster <- res_group$cluster_info %>%
    left_join(df_inst_counts) %>% arrange(cluster)
  ids.final <- res_group$selected_traits
  if(extra_traits != "None"){
    ids.final <- res_cluster %>% filter(id %in% extra_traits) %>% 
      bind_rows(res_cluster %>% filter(id %in% ids.final)) %>%
      distinct(cluster, .keep_all = T) %>% arrange(cluster) %>% pull(id)
  }
  filter.trait <- id.list[!id.list %in% ids.final]
  df_info <- left_join(df_info,res_cluster,by = c("id" = "id")) %>%
    rename(string_cluster = cluster)
  df_info[df_info$id %in% filter.trait,"status"] <- "delete due to string similarity"
  df_info[df_info$id %in% ids.final,"status"] <- "select after string similarity filtering"
  return(list(id.list=ids.final,trait.info=df_info))
}

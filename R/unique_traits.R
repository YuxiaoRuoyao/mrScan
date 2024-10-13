#' @title Calculate pairwise correlation and select unique traits
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param R_matrix Pairwise genetic correlation matrix. Colnames are trait names
#' @param df_pairs A dataframe contain string correlation for each pair
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @param method Filtering duplicate method: sample_size, nsnp, cluster. Default = "cluster"
#' @param cluster_selection_method Trait selection method in each cluster: n_inst (select traits with
#' the largest number of instruments), pvalue (-log10(p) mean for X and Y). Default = "n_inst"
#' @param extra_traits trait ID which is adjusted for in bidirection MR. Default = "None"
#' @param df_bidirection Dataframe of bidirection result. It should be input if you choose pvalue method for cluster selection.
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom dplyr left_join group_by slice_max pull arrange ungroup select
#' @export
unique_traits <- function(id.list,df_info,R_matrix,df_pairs,R2_cutoff=0.9,
                          method = "cluster",cluster_selection_method = "n_inst",
                          extra_traits = "None",df_bidirection = NULL){
  if(method=="cluster"){
    clusters <- greedy_cluster(id.list = id.list,R = R_matrix,
                               R2_cutoff = R2_cutoff,df_info = df_info)
    df_info <- left_join(df_info,clusters,by = c("id" = "id"))
    if(cluster_selection_method == "n_inst"){
      df_cluster <- left_join(clusters,df_info[,c("id","n_inst")],
                                     by=c("id" = "id")) %>% arrange(cluster)
      df_select <- df_cluster %>% group_by(cluster) %>%
        slice_max(n_inst, with_ties = FALSE) %>% ungroup()
    }else if(cluster_selection_method == "pvalue"){
      if(is.null(df_bidirection)){
        stop("Please input bidirection dataframe for pvalue method!")
      }else{
        df_cluster <- left_join(clusters,df_bidirection,by = c("id" = "id")) %>%
          mutate(mean_logp = rowMeans(cbind(-log10(p_ZtoY_adj), -log10(p_ZtoX_adj)))) %>%
          arrange(cluster)
        df_select <- df_cluster %>% group_by(cluster) %>% slice_max(mean_logp,with_ties = FALSE) %>%
          ungroup()
      }
    }
    if(extra_traits != "None"){
      ids.final <- df_cluster %>% filter(id %in% extra_traits) %>% bind_rows(df_select) %>%
        distinct(cluster, .keep_all = T) %>% arrange(cluster) %>% pull(id)
    }else{
      ids.final <- df_select %>% pull(id)
    }
    filter.trait <- df_cluster$id[!df_cluster$id %in% ids.final]
  }
  if(method=="nsnp"){
    df_pairs_snp <- df_pairs %>%
      left_join(df_info[,c("id","nsnp")],by=c("id1" = "id")) %>%
      left_join(df_info[,c("id","nsnp")],by=c("id2" = "id"),suffix = c("_id1", "_id2")) %>%
      left_join(df_info[,c("id","sample_size")],by=c("id1" = "id")) %>%
      left_join(df_info[,c("id","sample_size")],by=c("id2" = "id"),suffix = c("_id1", "_id2"))
    df_pairs_snp <- df_pairs_snp %>%
      group_by(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_")) %>%
      slice(1) %>% ungroup() %>% select(-grp)
    df_filter_snp <- df_pairs_snp[which(df_pairs_snp$cor>R2_cutoff),]
    filter.trait1 <- df_filter_snp %>% filter(nsnp_id1 > nsnp_id2) %>% pull(id2)
    filter.trait2 <- df_filter_snp %>% filter(nsnp_id1 < nsnp_id2) %>% pull(id1)
    filter.trait3 <- df_filter_snp %>% filter(nsnp_id1 == nsnp_id2) %>%
      filter(sample_size_id1 >= sample_size_id2) %>% pull(id2)
    filter.trait4 <- df_filter_snp %>% filter(nsnp_id1 == nsnp_id2) %>%
      filter(sample_size_id1 < sample_size_id2) %>% pull(id1)
    # If two traits have equal number of instruments, then select one with higher sample size
    filter.trait<-unique(c(filter.trait1,filter.trait2,filter.trait3,filter.trait4))
    ids.final <- id.list[!id.list %in% filter.trait]
  }
  if(method=="sample_size"){
    df_pairs_ss <- df_pairs %>%
      left_join(df_info[,c("id","nsnp")],by=c("id1" = "id")) %>%
      left_join(df_info[,c("id","nsnp")],by=c("id2" = "id"),suffix = c("_id1", "_id2")) %>%
      left_join(df_info[,c("id","sample_size")],by=c("id1" = "id")) %>%
      left_join(df_info[,c("id","sample_size")],by=c("id2" = "id"),suffix = c("_id1", "_id2"))
    df_pairs_ss <- df_pairs_ss %>%
      group_by(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_")) %>%
      slice(1) %>% ungroup() %>% select(-grp)
    df_filter_ss <- df_pairs_ss[which(df_pairs_ss$cor>R2_cutoff),]
    filter.trait1 <- df_filter_ss %>% filter(sample_size_id1 > sample_size_id2) %>% pull(id2)
    filter.trait2 <- df_filter_ss %>% filter(sample_size_id1 < sample_size_id2) %>% pull(id1)
    filter.trait3 <- df_filter_ss %>% filter(sample_size_id1 == sample_size_id2) %>%
      filter(nsnp_id1 >= nsnp_id2) %>% pull(id2)
    filter.trait4 <- df_filter_snp %>% filter(sample_size_id1 == sample_size_id2) %>%
      filter(nsnp_id1 < nsnp_id2) %>% pull(id1)
    # If two traits have equal number of instruments, then select one with higher sample size
    filter.trait<-unique(c(filter.trait1,filter.trait2,filter.trait3,filter.trait4))
    ids.final <- id.list[!id.list %in% filter.trait]
  }
  df_info[df_info$id %in% filter.trait,"status"] <- "delete due to duplicates"
  df_info[df_info$id %in% ids.final,"status"] <- "select after filtering duplicates"
  return(list(id.list=ids.final,trait.info=df_info))
}

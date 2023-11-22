#' @title Calculate pairwise correlation and select unique traits
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param method filtering duplicate method: sample_size, nsnp, cluster
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @importFrom dplyr left_join group_by pull slice_max
#' @export
unique_traits <- function(id.list,df_info,R2_cutoff=0.9,method){
  df_pairs<-calculate_cor(ids1 = id.list, ids2 = id.list)
  if(method=="cluster"){
    df_matrix <- as.data.frame.matrix(xtabs(cor ~ ., df_pairs))
    df_matrix <- abs(df_matrix)
    clusters <- greedy_cluster(id.list = names(df_matrix),R = df_matrix,
                               R2_cutoff = R2_cutoff)
    df_info <- dplyr::left_join(df_info,clusters,by = c("id" = "id"))
    df_cluster <- dplyr::left_join(clusters,df_info[,c("id","sample_size")],
                            by=c("id" = "id")) %>% dplyr::arrange(cluster)
    ids.final <- df_cluster %>% dplyr::group_by(cluster) %>%
      dplyr::slice_max(sample_size, with_ties = FALSE) %>% dplyr::ungroup() %>%
      dplyr::pull(id)
    filter.trait <- df_cluster$id[!df_cluster$id %in% ids.final]
  }
  if(method=="nsnp"){
    df_pairs_snp <- df_pairs %>%
      dplyr::left_join(df_info[,c("id","nsnp")],by=c("id1" = "id")) %>%
      dplyr::left_join(df_info[,c("id","nsnp")],by=c("id2" = "id"),suffix = c("_id1", "_id2")) %>%
      dplyr::left_join(df_info[,c("id","sample_size")],by=c("id1" = "id")) %>%
      dplyr::left_join(df_info[,c("id","sample_size")],by=c("id2" = "id"),suffix = c("_id1", "_id2"))
    df_pairs_snp <- df_pairs_snp %>%
      dplyr::group_by(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_")) %>%
      slice(1) %>% ungroup() %>% dplyr::select(-grp)
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
      dplyr::left_join(df_info[,c("id","nsnp")],by=c("id1" = "id")) %>%
      dplyr::left_join(df_info[,c("id","nsnp")],by=c("id2" = "id"),suffix = c("_id1", "_id2")) %>%
      dplyr::left_join(df_info[,c("id","sample_size")],by=c("id1" = "id")) %>%
      dplyr::left_join(df_info[,c("id","sample_size")],by=c("id2" = "id"),suffix = c("_id1", "_id2"))
    df_pairs_ss <- df_pairs_ss %>%
      dplyr::group_by(grp = paste(pmax(id1, id2), pmin(id1, id2), sep = "_")) %>%
      slice(1) %>% ungroup() %>% dplyr::select(-grp)
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

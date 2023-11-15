#' @title Calculate pairwise correlation and select unique traits
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param method filtering duplicate method: sample_size, nsnp, cluster
#' @param cor_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @param cluster_method the same with hclust method
#' @param ncluster The number of clusters. If not input by users, we'll use gap stats to
#' determine the optimal number of clusters
#' @param cluster_plot whether plot or not
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import mrScan
#' @import ieugwasr
#' @import dplyr
#' @import cluster
#' @import factoextra
#' @importFrom stats cutree dist hclust kmeans rect.hclust xtabs
#' @export
unique_traits <- function(id.list,df_info,cor_cutoff=0.9,method,cluster_method="complete",
                          ncluster=NULL,cluster_plot=FALSE){
  df_pairs<-calculate_cor(ids1 = id.list, ids2 = id.list)
  if(method=="cluster"){
    df_matrix <- as.data.frame.matrix(xtabs(cor ~ ., df_pairs))
    diag(df_matrix) <- 1
    df_matrix <- abs(df_matrix)
    clusters <- hclust(dist(df_matrix),method = cluster_method)
    plot(clusters)
    if(is.null(ncluster)){
      gap_stat <- cluster::clusGap(scale(df_matrix),
                                   FUN = kmeans, K.max = length(id.list)-1)
      # fviz_gap_stat(gap_stat)
      ncluster <- with(gap_stat,cluster::maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
      clusterCut <- data.frame(cutree(clusters, k=ncluster))
    }else{
      clusterCut <- data.frame(cutree(clusters, k=ncluster))
    }
    if(cluster_plot){
      plot(x = clusters, cex = 0.5)
      rect.hclust(tree = clusters, k = ncluster, which = 1:ncluster, border = 1:ncluster,
                      cluster = cutree(clusters, k = ncluster))
    }
    colnames(clusterCut)<-"cluster"
    clusterCut$id<-rownames(clusterCut)
    clusterCut <- clusterCut %>% left_join(df_info[,c("id","sample_size")],
                                           by=c("id" = "id")) %>% arrange(cluster)
    ids.final <- clusterCut %>% group_by(cluster) %>%
      slice_max(sample_size, with_ties = FALSE) %>% ungroup() %>% pull(id)
    filter.trait <- clusterCut$id[! clusterCut$id %in% ids.final]
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
    df_filter_snp <- df_pairs_snp[which(df_pairs_snp$cor>cor_cutoff),]
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
    df_filter_ss <- df_pairs_ss[which(df_pairs_ss$cor>cor_cutoff),]
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

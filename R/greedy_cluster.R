#' @title Greedy algorithem to get clusters
#' @param id.list GWAS ID list of traits based on previous steps
#' @param R Dataframe of pairwise correlation matrix. The colnames and rownames are trait IDs.
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @export
greedy_cluster <- function(id.list,R,R2_cutoff = 0.9){
  R <- R[id.list,id.list]
  groups <- list()
  unsorted <- seq(1,length(id.list))
  while(length(unsorted)>0){
    i <- length(groups)+1
    groups[[i]] <- unsorted[1]
    groupsize <- length(groups[[i]])
    done <- FALSE
    while (!done) {
      j <- apply(R[groups[[i]], ],2,max)
      k <- which(j > R2_cutoff)
      groups[[i]] <- unique(c(groups[[i]],k))
      groupsize_new <- length(groups[[i]])
      if(groupsize_new == groupsize){
        done <- TRUE
      }
      groupsize <- groupsize_new
    }
    unsorted <- unsorted[!unsorted %in% groups[[i]]]
  }
  cluster_info <- data.frame()
  for (m in 1:length(groups)) {
    sub <- id.list[groups[[m]]]
    df_sub <- data.frame(id = sub, cluster = m)
    cluster_info <- rbind(cluster_info,df_sub)
  }
  return(cluster_info)
}

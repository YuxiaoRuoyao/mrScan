#' @title Greedy algorithm to get clusters.
#' @description Hierarchically clusters traits using a greedy algorithm.
#' Initially, traits are clustered according to the highest pairwise correlations
#' above a cutoff threshold. Then, within each cluster, traits are subclustered
#' based on descending order of instrument numbers with high correlations grouped together.
#' @param id.list GWAS ID list of traits based on previous steps
#' @param R Dataframe of pairwise correlation matrix. The colnames and rownames are trait IDs.
#' @param df_info Dataframe of trait ID info. It requires to contain `id`: trait ID,
#' `n_inst`: the number of instruments.
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @export
greedy_cluster <- function(id.list, R, df_info, R2_cutoff = 0.9) {
  R <- R[id.list, id.list]
  groups <- list()
  unsorted <- seq(1, length(id.list))
  while (length(unsorted) > 0) {
    i <- length(groups) + 1
    groups[[i]] <- unsorted[1]
    groupsize <- length(groups[[i]])
    done <- FALSE
    while (!done) {
      j <- apply(R[groups[[i]], ], 2, max)
      k <- which(j > R2_cutoff)
      groups[[i]] <- unique(c(groups[[i]], k))
      groupsize_new <- length(groups[[i]])
      if (groupsize_new == groupsize) {
        done <- TRUE
      }
      groupsize <- groupsize_new
    }
    unsorted <- unsorted[!unsorted %in% groups[[i]]]
  }
  cluster_info <- data.frame()
  num_instruments <- setNames(df_info$n_inst[match(id.list, df_info$id)], id.list)
  for (i in seq_along(groups)) {
    subcluster <- 1
    remaining <- groups[[i]]
    while(length(remaining) > 0) {
      current_id_index <- remaining[which.max(num_instruments[id.list[remaining]])]
      current_id <- id.list[current_id_index]
      correlated <- which(R[current_id, id.list[remaining]] > R2_cutoff)
      correlated_ids <- id.list[remaining[correlated]]
      cluster_info <- rbind(cluster_info, data.frame(
        id = correlated_ids,
        cluster = i,
        subcluster = subcluster
      ))
      remaining <- remaining[-correlated]
      subcluster <- subcluster + 1
    }
  }
  return(cluster_info)
}

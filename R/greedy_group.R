#' @title Greedy algorithm to get unique traits.
#' @description Hierarchically clusters traits using a greedy algorithm.
#' First rank traits according to the number of genome-wide significant variants after LD clumping. 
#' We then iteratively select the top trait and remove all traits with a similarity score greater than 
#' the threshold with the selected trait, continuing until all traits are either selected or removed.
#' @param id.list GWAS ID list of traits based on previous steps
#' @param R Dataframe of pairwise correlation matrix. The colnames and rownames are trait IDs.
#' @param df_info Dataframe of trait ID info. It requires to contain `id`: trait ID,
#' `n_inst`: the number of instruments.
#' @param R2_cutoff high correlation cutoff to assign as duplicated traits. Default=.9
#' @export
greedy_group <- function(id.list, R, df_info, R2_cutoff = 0.9) {
  R <- R[id.list, id.list, drop = FALSE]
  num_instruments <- setNames(df_info$n_inst[match(id.list, df_info$id)], id.list)
  if (any(is.na(num_instruments))) {
    stop("Some traits in id.list are missing n_inst in df_info.")
  }
  cluster_info <- data.frame()
  selected_traits <- character(0)
  remaining <- id.list
  group_idx <- 1
  while (length(remaining) > 0) {
    # select the remaining trait with the largest n_inst
    current_id <- remaining[which.max(num_instruments[remaining])]
    selected_traits <- c(selected_traits, current_id)
    # group all remaining traits above the threshold with the selected trait
    grouped_ids <- remaining[which(R[current_id, remaining] > R2_cutoff)]
    cluster_info <- rbind(
      cluster_info,
      data.frame(
        id = grouped_ids,
        cluster = group_idx,
        stringsAsFactors = FALSE))
    # remove grouped traits and continue
    remaining <- remaining[!remaining %in% grouped_ids]
    group_idx <- group_idx + 1
  }
  return(list(cluster_info = cluster_info, selected_traits = selected_traits))
}
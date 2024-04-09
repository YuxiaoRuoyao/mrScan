#' @title Format and combine GWAS summary data, do LDSC for R matrix and LD-pruning
#' @description This is a combined function for a series of step: 1. Format and combine
#' GWAS summary statistics; 2. Do LD pruning for every chromosome and combine them
#' together; 3. Use LDSC to calculate sample overlap and genetic correlation matrix.
#' We recommend you to use this only when the trait list is short.
#' @param df_file A dataframe contain "id" (trait id) and "location" (GWAS summary data location) columns.
#' @param df_info Dataframe of trait info from previous step
#' @param r2_thresh Clumping r2 cut off. Default = 0.001
#' @param clump_kb Clumping distance cut off. Default = 10000
#' @param type LD clumping prioritization. Either pvalue or rank. Default = "pvalue"
#' @param pthresh pvalue threshold. Default = 1
#' @param ref_path Path for the LD reference panel.
#' @param ld_files Paths of reference LD score files
#' @param m_files Paths of reference M files
#' @returns A list of selected traits, a dataframe of conditional instrument strength and a dataframe of trait info
#'
#' @import dplyr
#' @import MVMR
#' @importFrom purrr map_dfr
#' @export
format_ldsc_prune <- function(df_file,df_info,r2_thresh,clump_kb,type,
                              pthresh,ref_path,ld_files,m_files){
  formatted_dats <- list()
  pruned_dats <- list()
  for (c in 1:22) {
    formatted_dat <- format_combine_gwas(df_file = df_file, c = c, df_info = df_info)
    pruned_dat <- ld_prune_plink(X = formatted_dat, r2_thresh = r2_thresh, clump_kb = clump_kb,
                                 ref_path = ref_path, type = type, pthresh = pthresh)
    pruned_dats[[c]] <- pruned_dat
    formatted_dats[[c]] <- formatted_dat
  }
  combine_formatted <- do.call(rbind, formatted_dats)
  R <- ldsc_full(dat = combine_formatted, ld_files = ld_files, m_files = m_files)
  combine_pruned <- do.call(rbind, pruned_dats)
  return(list(R = R, dat = combine_pruned))
}

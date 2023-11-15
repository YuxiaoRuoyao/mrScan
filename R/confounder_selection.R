#' @title Select confounders by variable selection
#' @param id_exposure GWAS ID of the main exposure
#' @param id_outcome GWAS ID of the outcome
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param method Bidirection MR method. Same with TwoSampleMR package.Default = "mr_raps"
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import mrScan
#' @import ieugwasr
#' @import dplyr
#' @import data.table
#' @importFrom stats pnorm
#' @export
confounder_selection <- function(){

}

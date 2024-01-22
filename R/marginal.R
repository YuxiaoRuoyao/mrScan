#' @title Use marginal selection to do confounder selection
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @param df_bidirection Dataframe of trait with four direction estimates. Result from downstream_filter
#' @param extra_traits trait ID which is adjusted for in bidirection MR. Default = "None"
#' @param p_cutoff pvalue threshold. Default = 0.05
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import dplyr
#' @export
marginal <- function(id.list,df_info,df_bidirection,extra_traits = "None",
                     p_cutoff = 0.05){
  df_bidirection <- df_bidirection %>%
    mutate(p_ZtoX = 2*(1-pnorm(abs(b_ZtoX/se_ZtoX))),
           p_ZtoY = 2*(1-pnorm(abs(b_ZtoY/se_ZtoY))))
  id.select <- df_bidirection %>% filter(id %in% id.list) %>%
    filter(p_ZtoX < p_cutoff & p_ZtoY < p_cutoff) %>% pull(id)
  if(extra_traits != "None"){
    id.select <- c(id.select,extra_traits)
  }
  df_info[df_info$id %in% id.select,"status"] <- "Select by marginal selection"
  return(list(id.list=id.select,trait.info=df_info))
}

#' @title Calculate trait similarity by trait name
#' @param id.list GWAS ID list of traits based on previous steps
#' @param df_info Dataframe of trait info from previous steps
#' @returns df_pair: A dataframe contain string correlation for each pair. R_matrix: Pairwise correlation matrix
#'
#' @import RecordLinkage
#' @import dplyr
#' @importFrom purrr map2 map
#' @export
#'
string_sim <- function(id.list,df_info){
  df_pairs <- data.frame(t(combn(id.list, 2)))
  df_trait <- df_info[,c("id","trait")] %>% filter(id %in% id.list)
  A<-strsplit(df_trait$trait, '[()]')
  df_trait$item1 <- unlist(map(A,1))
  item2 <- map(A,2)
  item2[sapply(item2, is.null)] <- NA
  df_trait$item2 <- unlist(item2)
  df_pairs <- left_join(df_pairs,df_trait,by = c("X1" = "id")) %>%
    rename("trait1" = "trait","trait1_item1" = "item1", "trait1_item2" = "item2") %>%
    left_join(df_trait,by = c("X2" = "id")) %>%
    rename("trait2" = "trait","trait2_item1" = "item1","trait2_item2" = "item2")
  df_pairs <- df_pairs %>% mutate_all(tolower)
  df_pairs$string_cor1 <- map2(df_pairs$trait1_item1, df_pairs$trait2_item1, function(x, y){
    jarowinkler(x,y)
  }) %>% unlist()
  df_pairs$string_cor2 <- map2(df_pairs$trait1_item1, df_pairs$trait2_item2, function(x, y){
    jarowinkler(x,y)
  }) %>% unlist()
  df_pairs$string_cor3 <- map2(df_pairs$trait2_item1, df_pairs$trait1_item2, function(x, y){
    jarowinkler(x,y)
  }) %>% unlist()
  df_pairs$string_cor <- do.call(pmax, c(df_pairs[,c("string_cor1","string_cor2","string_cor3")],
                                         list(na.rm=TRUE)))
  df_pairs <- df_pairs %>% dplyr::select(X1,X2,trait1,trait2,string_cor)
  res <- df_pairs %>% dplyr::select(X1,X2,string_cor)
  res2 <- data.frame(X1 = res$X2, X2 = res$X1, string_cor = res$string_cor)
  res_all <- rbind(res,res2)
  df_matrix <- as.data.frame.matrix(xtabs(string_cor ~ ., res_all))
  df_matrix <- abs(df_matrix)
  return(list(df_pair = df_pairs, R_matrix = df_matrix))
}

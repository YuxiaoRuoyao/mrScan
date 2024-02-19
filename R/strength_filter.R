#' @title Calculates the conditional F-statistic to assess instrument strength and filter traits
#' @param beta_files Paths of merged GWAS summary data after LD pruning
#' @param R_matrix Pairwise sample overlap matrix among traits
#' @param df_info Dataframe of trait info from previous step
#' @param pval_threshold pvalue cutoff for selecting instruments. Default = 5e-8
#' @param F_threshold F-statistic cutoff. Default = 5
#' @param extra_traits trait ID you want to include no matter the instrument strength. Default = "None"
#' @returns A list of selected traits, a dataframe of conditional instrument strength and a dataframe of trait info
#'
#' @import dplyr
#' @import MVMR
#' @importFrom purrr map_dfr
#' @export
strength_filter <- function(beta_files,R_matrix,df_info,
                            pval_threshold = 5e-8,F_threshold = 5, extra_traits = "None"){
  X <- purrr::map_dfr(beta_files, readRDS)
  snp <- data.frame(X$snp)
  beta_hat <- X %>% select(ends_with(".beta"))
  se <- X %>% select(ends_with(".se"))
  p <- X %>% select(ends_with(".p"))
  nms <- stringr::str_replace(names(beta_hat), ".beta", "")
  names(beta_hat)<-names(se)<-names(p)<-nms
  o <- match(colnames(R_matrix), nms)
  beta_hat <- data.frame(beta_hat[, o],check.names = F)
  se <- data.frame(se[, o],check.names = F)
  i <- ncol(beta_hat)
  pmin <- apply(p[,-1, drop = F], 1, min)
  ix <- which(pmin < pval_threshold)
  exp_se <- as.matrix(se[ix,2:i])
  F.data<- format_mvmr(BXGs = as.matrix(beta_hat[ix, 2:i]),
                       BYG = beta_hat[ix,1],
                       seBXGs = exp_se,
                       seBYG = se[ix,1],
                       RSID = snp[ix,1])
  sigmalist <- vector("list", length(ix))
  omega <- R_matrix[2:i,2:i]
  for (j in 1:length(ix)) {
    se_matrix <- diag(exp_se[j,],nrow = i-1)
    sigmalist[[j]] <- se_matrix %*% omega %*% se_matrix
  }
  sres <- data.frame(t(strength_mvmr(r_input = F.data, gencov = sigmalist)))
  sres$id<-colnames(beta_hat[-1])
  select.id <- sres %>% filter(F.statistic > F_threshold) %>% pull(id)
  if(extra_traits != "None"){
    select.id <- c(select.id,extra_traits)
  }
  other.id <- sres$id[!sres$id %in% select.id]
  df_info[df_info$id %in% other.id,"status"] <- "Delete due to weak instruments strength"
  return(id.list = select.id,df_strength = data.frame(sres),trait.info=df_info)
}

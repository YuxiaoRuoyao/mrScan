#' @title Estimate genetic correlation matrix by pvalue
#' @param beta_dir Directory path of merged data. Default is the current work directory.
#' @param ref_path Path for the LD reference panel.
#' @param out_dir Output data path. Default is in the current work directory.
#' @param prefix Name prefix for the output. Default = NULL
#' @returns Save one dataframe per chromosome with columns for SNP info
#'
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import sumstatFactors
#' @export
estimate_R_pval <- function(ld_prune_file_dir=NULL,prefix=NULL,p_thresh=0.05,
                            out_dir = NULL){
  files <- paste0(ld_prune_file_dir,prefix,".beta.ldpruned.",seq(1,22),".RDS")
  X <- purrr::map_dfr(files, readRDS)
  beta_hat <- X %>%
    select(ends_with(".beta")) %>%
    as.matrix()
  se_hat <- X %>%
    select(ends_with(".se")) %>%
    as.matrix()
  nms <- colnames(beta_hat) %>% stringr::str_replace(".beta$", "")
  Rpt <- R_pt(B_hat = beta_hat,
              S_hat = se_hat,
              p_val_thresh = p_thresh,
              return_cor = TRUE,
              make_well_conditioned = FALSE
  )
  colnames(Rpt) <- rownames(Rpt) <- nms
  saveRDS(Rpt, file=paste0(out_dir,prefix,".R_est_pval.RDS"))
}

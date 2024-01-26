#' @title Initially extract candidate traits
#' @param id_exposure GWAS ID for main exposure.
#' @param pval_x p-value threshold to extract instruments of the main exposure. Default=5e-8
#' @param pval_z p-value threshold to retrieve traits. Default=1e-5
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are "EUR","SAS","EAS","AFR","AMR"
#' @param batch GWAS database sub-batches, a vector. Default = c("ieu-a", "ieu-b","ukb-b")
#' @param r2 clumping r2 threshold. Default=0.001
#' @param kb clumping kb window. Default=10000
#' @param access_token Google OAuth2 access token. Default=check_access_token()
#' @param min_snps the number of minimum shared SNPs with IV of X. Default=5
#' @param type_exposure Exposure data type. Either could be "IEU" or "local". Default = "IEU"
#' @param type_candidate_traits Candidate traits data type. Either could be "IEU" or "local". Default = "IEU"
#' @param file_path File path of local exposure GWAS summary data. It should be entered when you use local exposure data. Default = NULL
#' @param ref_path LD reference data path. It should be entered when you use local exposure data. Default = NULL
#' @param file_list GWAS summary data file paths for candidate traits. It should be entered when you use local traits data. Default = NULL
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import TwoSampleMR
#' @import ieugwasr
#' @import dplyr
#' @import stringr
#' @import readr
#' @import gwasvcf
#' @import rlang
#' @export
retrieve_traits <- function (id_exposure, pval_x = 5e-8, pval_z = 1e-5,
                             pop = "EUR", batch = c("ieu-a", "ieu-b","ukb-b"),
                             r2 = 0.001, kb = 10000,
                             access_token = ieugwasr::check_access_token(),
                             min_snps = 5,
                             type_exposure = "IEU",
                             type_candidate_traits = "IEU",
                             file_path = NULL, ref_path = NULL,file_list = NULL,
                             trait_list = NULL, snp_name_list = NULL,
                             beta_hat_name_list = NULL, se_name_list = NULL,
                             p_value_name_list = NULL) {
  df_inst <- get_exposure_inst(id_x = id_exposure,type = type_exposure,
                               file_path = file_path,
                               pval_x = 5e-8,r2 = r2,kb = kb, pop = pop,
                               access_token = access_token,ref_path = ref_path)
  df_association <- get_association_inst(df_inst = df_inst,type = type_candidate_traits,
                                         pval_z = pval_z,batch = batch,
                                         access_token = access_token,file_list = file_list,
                                         trait_list = trait_list,snp_name_list= snp_name_list,
                                         beta_hat_name_list = beta_hat_name_list,
                                         se_name_list = se_name_list,
                                         p_value_name_list = p_value_name_list)
  x <- df_association %>% dplyr::group_by(id) %>% dplyr::summarize(n = length(unique(rsid)))
  cat(sum(x$n >= min_snps), "traits have at least", min_snps,
      "shared variants with", id_exposure, "\n")
  ids <- x$id[x$n >= min_snps]
  df_trait <- gwasinfo(ids)
  df_trait <- df_trait[,c("id","trait","sex","consortium","nsnp","note","sample_size",
                          "population","year")]
  df_trait['status'] <- 'Initial List'
  # Delete X
  id.list.initial <- ids[!ids %in% id_exposure]
  df_trait[df_trait$id %in% id_exposure,"status"] <- "Delete due to it's the main exposure"
  return(list(id.list=id.list.initial,trait.info=df_trait))
}

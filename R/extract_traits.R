#' @title Initial extract and filtering
#' @param id_exposure GWAS ID for main exposure.
#' @param id_outcome GWAS ID for the outcome.
#' @param batch GWAS database sub-batches, a vector.
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are "EUR","SAS","EAS","AFR","AMR"
#' @param pval_x p-value threshold to extract instruments of the main exposure, default=5e-8
#' @param pval_z p-value threshold to retreive traits, default=1e-5
#' @param r2 clumping r2 threshold. Default=0.001
#' @param kb clumping kb window. Default=10000
#' @param access_token Google OAuth2 access token. Default=check_access_token()
#' @param min_snps the number of minimum shared SNPs with IV of X. Default=5
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import mrScan
#' @import ieugwasr
#' @import dplyr
#' @export
extract_traits<-function(id_exposure,id_outcome,batch = c("ieu-a", "ieu-b","ukb-b","ebi-a"),
                         pop = "EUR",pval_x=5e-8, pval_z=1e-5,r2=0.001,kb = 10000,
                         access_token = ieugwasr::check_access_token(),min_snps=5){
  phe <- mrScan::retrieve_traits(id_exposure, pval_x,pval_z,pop = pop, batch=batch,
                                 r2 = r2, kb = kb,
                                 access_token = access_token,min_snps =min_snps)
  id.list <- unique(phe$phe$id)
  # Delete X
  id.list.initial <- id.list[!id.list %in% id_exposure]
  df_trait <- ieugwasr::gwasinfo(id.list.initial)
  df_trait <- df_trait[,c("id","trait","sex","consortium","nsnp","note","sample_size",
                          "pmid","population","year")]
  df_trait['status'] <- 'Initial List'
  return(list(id.list=id.list.initial,trait.info=df_trait))
}

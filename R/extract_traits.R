#' @title Initial extract and filtering
#' @param id_exposure GWAS ID for main exposure.
#' @param id_outcome GWAS ID for the outcome.
#' @param batch GWAS database sub-batches, a vector.
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are "EUR","SAS","EAS","AFR","AMR"
#' @param pval_x p-value threshold to extract instruments of the main exposure, default=5e-8
#' @param pval_z p-value threshold to retrieve traits, default=1e-5
#' @param r2 clumping r2 threshold. Default=0.001
#' @param kb clumping kb window. Default=10000
#' @param access_token Google OAuth2 access token. Default=check_access_token()
#' @param min_snps the number of minimum shared SNPs with IV of X. Default=5
#' @param min_instruments the number of minimum instruments for each selected traits. Default=3
#' @returns A GWAS ID vector and a trait info dataframe
#'
#' @import ieugwasr
#' @import dplyr
#' @export
extract_traits<-function(id_exposure,id_outcome,batch = c("ieu-a","ieu-b","ukb-b"),
                         pop = "EUR",pval_x=5e-8, pval_z=1e-5,r2=0.001,kb = 10000,
                         access_token = ieugwasr::check_access_token(),min_snps=5,
                         min_instruments = 3){
  batch1 <- c("ieu-a","ieu-b","ukb-b")
  batch2 <- batch[!batch %in% batch1]
  phe1 <- retrieve_traits(id_exposure, pval_x,pval_z,pop = pop, batch=batch1,
                          r2 = r2, kb = kb,
                          access_token = access_token,min_snps =min_snps)
  if(length(batch2) != 0){
    phe2 <- retrieve_traits(id_exposure, pval_x,pval_z,pop = pop, batch=batch2,
                            r2 = r2, kb = kb,
                            access_token = ieugwasr::check_access_token(),
                            min_snps =min_snps)
    id.list <- unique(c(phe1$phe$id,phe2$phe$id))
  }else{
    id.list <- unique(phe1$phe$id)
  }
  # Delete X
  id.list.initial <- id.list[!id.list %in% id_exposure]
  df_trait <- gwasinfo(id.list.initial)
  df_trait <- df_trait[,c("id","trait","sex","consortium","nsnp","note","sample_size",
                          "population","year")]
  df_trait['status'] <- 'Initial List'
  return(list(id.list=id.list.initial,trait.info=df_trait))
}

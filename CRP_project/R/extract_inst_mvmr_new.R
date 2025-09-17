library(TwoSampleMR)
library(ieugwasr)

x1 <- snakemake@params[["trait1"]]
x2 <- snakemake@params[["trait2"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out <- snakemake@output[["out"]]

df_inst <- mv_extract_exposures_local(
  filenames_exposure = c("/nfs/turbo/sph-jvmorr/CRP_project/GWAS_summary_data/GCST90475667.tsv.gz",
                         "/nfs/turbo/sph-jvmorr/CRP_project/GWAS_summary_data/ukb-b-19953.tsv.gz"),
  sep = "\t",
  snp_col = c("hm_rsid","rsid"),
  beta_col = c("hm_beta","ES"),
  se_col = c("standard_error","SE"),
  eaf_col = c("hm_effect_allele_frequency","AF"),
  effect_allele_col = c("hm_effect_allele","ALT"),
  other_allele_col = c("hm_other_allele","REF"),
  pval_col = "p_value",
  ncase_col = c("num_cases","ncase"),
  ncontrol_col = c("num_controls","ncontrol"),
  samplesize_col = c("sample_size","SS"),
  pval_threshold = inst_pval,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR"
)
df_inst$id.exposure[which(df_inst$exposure == "exposure1")] <- x1
df_inst$id.exposure[which(df_inst$exposure == "exposure2")] <- x2

saveRDS(df_inst,file = out)

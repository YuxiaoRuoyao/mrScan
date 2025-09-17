library(TwoSampleMR)
library(ieugwasr)
x <- snakemake@params[["trait"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
filename <- snakemake@params[["filename"]]
out <- snakemake@output[["out"]]

df_inst <- read_exposure_data(
    filename = filename,
    sep = "\t",
    snp_col = "hm_rsid",
    beta_col = "hm_beta",
    se_col = "standard_error",
    effect_allele_col = "hm_effect_allele",
    other_allele_col = "hm_other_allele",
    eaf_col = "hm_effect_allele_frequency",
    pval_col = "p_value",
    ncase_col = "num_cases", 
    ncontrol_col = "num_controls",
    samplesize_col = "sample_size",
    chr_col = "hm_chrom", pos_col = "hm_pos"
)
df_inst <- clump_data(df_inst,plink_bin = genetics.binaRies::get_plink_binary(), 
                      bfile = "/nfs/turbo/sph-jvmorr/ld_reference_files/1kg_plink/EUR")
df_inst$exposure <- "T2D"
df_inst$id.exposure <- x
df_inst <- df_inst[df_inst$pval.exposure < inst_pval,]
saveRDS(df_inst,file = out)

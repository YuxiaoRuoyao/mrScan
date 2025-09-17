library(TwoSampleMR)
library(ieugwasr)

x1 <- snakemake@params[["trait1"]]
x2 <- snakemake@params[["trait2"]]
x3 <- snakemake@params[["trait3"]]
id_outcome <- snakemake@params[["id_outcome"]]
filename <- snakemake@params[["filename"]]
inst_pval <- as.numeric(snakemake@params[["pval_instruments"]])
out1 <- snakemake@output[["out1"]]
out2 <- snakemake@output[["out2"]]
out3 <- snakemake@output[["out3"]]
mv_extract_exposures_edit <- function (id_exposure, clump_r2 = 0.001, clump_kb = 10000, harmonise_strictness = 2,
          opengwas_jwt = ieugwasr::get_opengwas_jwt(), find_proxies = TRUE,
          force_server = FALSE, pval_threshold = 5e-08, pop = "EUR",
          plink_bin = NULL, bfile = NULL)
{
  stopifnot(length(id_exposure) > 1)
  id_exposure <- ieugwasr::legacy_ids(id_exposure)
  exposure_dat <- extract_instruments(id_exposure, p1 = pval_threshold,
                                      r2 = clump_r2, kb = clump_kb, opengwas_jwt = opengwas_jwt,
                                      force_server = force_server)
  temp <- exposure_dat
  temp$id.exposure <- 1
  temp <- temp[order(temp$pval.exposure, decreasing = FALSE),
  ]
  temp <- subset(temp, !duplicated(SNP))
  temp <- clump_data(temp, clump_p1 = pval_threshold, clump_r2 = clump_r2,
                     clump_kb = clump_kb, pop = pop, plink_bin = plink_bin,
                     bfile = bfile)
  exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)
  d1 <- extract_outcome_data(exposure_dat$SNP, id_exposure,
                             opengwas_jwt = opengwas_jwt, proxies = find_proxies, splitsize=50, proxy_splitsize=20)
  stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
  d1 <- subset(d1, mr_keep.outcome)
  d2 <- subset(d1, id.outcome != id_exposure[1])
  d1 <- convert_outcome_to_exposure(subset(d1, id.outcome ==
                                             id_exposure[1]))
  d <- harmonise_data(d1, d2, action = harmonise_strictness)
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure) - 1]
  d <- subset(d, SNP %in% keepsnps)
  dh1 <- subset(d, id.outcome == id.outcome[1], select = c(SNP,
                                                           exposure, id.exposure, effect_allele.exposure, other_allele.exposure,
                                                           eaf.exposure, beta.exposure, se.exposure, pval.exposure))
  dh2 <- subset(d, select = c(SNP, outcome, id.outcome, effect_allele.outcome,
                              other_allele.outcome, eaf.outcome, beta.outcome, se.outcome,
                              pval.outcome))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  dh <- rbind(dh1, dh2)
  return(dh)
}
ex_dat <- mv_extract_exposures_edit(id_exposure = c(x1,x2,x3), pval_threshold = inst_pval)
df <- data.table::fread(filename, header = TRUE)
#df$beta_hat <- log(df$odds_ratio)
#df$standard_error <- (log(df$ci_upper) - log(df$ci_lower)) / (2 * qnorm(0.975))
#df$p_value <- 2 * (1 - pnorm(abs(df$beta_hat/df$standard_error)))
#colnames(df)[1]="hm_chrom"
#colnames(df)[2]="hm_pos"
#colnames(df)[3]="hm_effect_allele"
#colnames(df)[4]="hm_other_allele"
#colnames(df)[7]="hm_effect_allele_frequency"
#colnames(df)[9]="hm_rsid"
#colnames(df)[13]<-"sample_size"
#colnames(df)[25]<-"hm_beta"
#subdf <- df %>% select(hm_rsid,hm_pos,hm_chrom,hm_effect_allele,hm_other_allele,
#                       hm_beta,standard_error,p_value,hm_effect_allele_frequency,
#                       sample_size,num_cases,num_controls)
#write.table(subdf,file=gzfile("/nfs/turbo/sph-jvmorr/CRP_project/GWAS_summary_data/GCST90475667.tsv.gz"),
#            sep = "\t", row.names = F)

out_dat <- format_data(as.data.frame(df), type = "outcome",
                       snps = ex_dat$SNP, snp_col = "hm_rsid",
                       beta_col = "hm_beta", se_col = "standard_error", eaf_col = "hm_effect_allele_frequency",
                       effect_allele_col = "hm_effect_allele", other_allele_col = "hm_other_allele",
                       pval_col = "p_value", ncase_col = "num_cases",
                       ncontrol_col = "num_controls", samplesize_col = "sample_size",
                       chr_col = "hm_chrom", pos_col = "hm_pos")
out_dat$data_source.outcome <- "textfile"
dat <- mv_harmonise_data(ex_dat,out_dat)
saveRDS(dat,file = out1)
saveRDS(out_dat,file = out2)
saveRDS(ex_dat,file = out3)

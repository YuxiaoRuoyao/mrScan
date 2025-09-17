library(dplyr)
library(mrScan)
library(MVMR)
library(ieugwasr)

dat <- snakemake@input[["dat"]]
#out_dat <- snakemake@input[["out_dat"]]
ex_dat <- snakemake@input[["ex_dat"]]
df_info <- read.csv(snakemake@input[["trait_info"]])
pval_threshold <- as.numeric(snakemake@params[["pval_threshold"]])
F_threshold <- as.numeric(snakemake@params[["F_threshold"]])
effect_size_cutoff <- as.numeric(snakemake@params[["effect_size_cutoff"]])
type_outcome <- snakemake@params[["type_outcome"]]
prevalence_outcome <- as.numeric(snakemake@params[["prevalence_outcome"]])
id_exposure <- snakemake@params[["id_exposure"]]
extra_trait <- snakemake@params[["extra_trait"]]
out <- snakemake@output[["out"]]

dat <- readRDS(dat)
ex_dat <- readRDS(ex_dat)
#out_dat <- readRDS(out_dat)
df_af_exp <- generate_df_af_exp(ex_dat = ex_dat,mv_dat = dat)
#df_af_out <- data.frame(SNP = rownames(dat$exposure_beta), beta = dat$outcome_beta) %>%
#  left_join(out_dat, by = "SNP") %>%
#  mutate(eaf.outcome = ifelse(abs(beta - beta.outcome) < 1e-8, eaf.outcome, 1 - eaf.outcome))

#res <- strength_filter(
#  dat = dat,
#  dat_type = "IEU",
#  df_info = df_info,
#  pval_threshold = pval_threshold,
#  F_threshold = F_threshold,
#  effect_size_cutoff = effect_size_cutoff,
#  Filter = FALSE,
#  type_outcome = type_outcome,
#  prevalence_outcome = prevalence_outcome,
#  df_af_exp = df_af_exp, df_af_out = df_af_out)
res <- strength_filter(
  dat = dat,
  dat_type = "IEU",
  df_info = df_info,
  pval_threshold = pval_threshold,
  F_threshold = F_threshold,
  effect_size_cutoff = effect_size_cutoff,
  Filter = FALSE,
  type_outcome = type_outcome,
  prevalence_outcome = prevalence_outcome, df_af_exp = df_af_exp)
if (extra_trait != "None") {
    variable_id_row <- which(!(res$id %in% c(id_exposure, extra_trait)))
    res <- res[c(variable_id_row, which(res$id == id_exposure), which(res$id == extra_trait)), ]
  } else {
    variable_id_row <- which(!(res$id %in% id_exposure))
    res <- res[c(variable_id_row, which(res$id == id_exposure)), ]
}
df <- data.frame(t(c(as.matrix(res))))
colnames(df) <- paste0(rep(colnames(res), each = nrow(res)), 1:nrow(res))

saveRDS(df,file = out)

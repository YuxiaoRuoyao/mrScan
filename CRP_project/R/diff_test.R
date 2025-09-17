library(dplyr)
outcomes <- list(
  list(id = "ebi-a-GCST005194", label = "CAD"),
  list(id = "ebi-a-GCST006906", label = "stroke"),
  list(id = "ebi-a-GCST90027158", label = "Alzheimer"),
  list(id = "ieu-b-7", label = "Parkinson"),
  list(id = "ieu-b-5102", label = "SCZ"),
  list(id = "ieu-b-41", label = "BD"),
  list(id = "ebi-a-GCST90038683", label = "IBD"),
  list(id = "ebi-a-GCST90018910", label = "RA"),
  list(id = "ebi-a-GCST90475667", label = "T2D"),
  list(id = "ebi-a-GCST90018808", label = "Colorectal"),
  list(id = "ebi-a-GCST007090", label = "Knee"),
  list(id = "ebi-a-GCST90014022", label = "BMD"),
  list(id = "ieu-b-109", label = "HDL"),
  list(id = "ieu-b-110", label = "LDL"),
  list(id = "ieu-b-111", label = "Triglycerides"),
  list(id = "ebi-a-GCST90014006", label = "HbA1c")
)
analyses <- list(
  list(suffix = "_MVMR_GRAPPLE_FDR_p_0.05_summary.csv",  analysis = "FDR_0.05"),
  list(suffix = "_MVMR_GRAPPLE_FDR_p_0.1_summary.csv",   analysis = "FDR_0.1"),
  list(suffix = "_MR_GRAPPLE_FDR_p_0.05_summary.csv",    analysis = "no_BMI")
)

df_all <- data.frame()
for (outcome in outcomes) {
  id <- outcome$id
  label <- outcome$label
  for (ana in analyses) {
    ana_path <- paste0(id, "/results/", label, ana$suffix)
    df_tmp <- read.csv(ana_path)
    df_tmp$outcome <- label
    df_tmp$analysis <- ana$analysis
    df_all <- bind_rows(df_all, df_tmp)
  }
}
df_all <- df_all %>% filter(exposure == "ebi-a-GCST90029070") %>%
    filter(method != "IVW_T_5e-08" & method != "ESMR_5e-08" & method != "ESMR_optimize_5e-08")
df_uvmr <- df_all %>% filter(type == "UVMR")
df_mvmr <- df_all %>% filter(type == "Stepwise")
df_joined <- inner_join(
  df_uvmr, df_mvmr,
  by = c("outcome", "method", "analysis"),
  suffix = c("_UVMR", "_MVMR")
) %>% 
    select(outcome, method, analysis, b_UVMR, b_MVMR, se_UVMR, se_MVMR) %>%
    mutate(
    Z_diff = (b_UVMR - b_MVMR) / sqrt(se_UVMR^2 + se_MVMR^2),
    p_diff = 2 * (1 - pnorm(abs(Z_diff)))
  )
write.csv(df_joined,"res_diff.csv",row.names=F)
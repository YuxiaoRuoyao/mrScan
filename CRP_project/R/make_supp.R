# Supplementary table making
library(openxlsx)
library(dplyr)
library(ieugwasr)

outcomes <- list(
  list(id = "ebi-a-GCST005194", label = "CAD", sheet_name = "Coronary Artery Disease"),
  list(id = "ebi-a-GCST006906", label = "stroke", sheet_name = "Stroke"),
  list(id = "ebi-a-GCST90027158", label = "Alzheimer", sheet_name = "Alzheimer's disease"),
  list(id = "ieu-b-7", label = "Parkinson", sheet_name = "Parkinson's disease"),
  list(id = "ieu-b-5102", label = "SCZ", sheet_name = "Schizophrenia"),
  list(id = "ieu-b-41", label = "BD", sheet_name = "Bipolar Disorder"),
  list(id = "ebi-a-GCST90038683", label = "IBD", sheet_name = "InflammatoryBowelDisease"),
  list(id = "ebi-a-GCST90018910", label = "RA", sheet_name = "Rheumatoid Arthritis"),
  list(id = "ebi-a-GCST90475667", label = "T2D", sheet_name = "Type2 Diabetes"),
  list(id = "ebi-a-GCST90018808", label = "Colorectal", sheet_name = "Colorectal Cancer"),
  list(id = "ebi-a-GCST007090", label = "Knee", sheet_name = "Knee Osteoarthritis"),
  list(id = "ebi-a-GCST90014022", label = "BMD", sheet_name = "Bone Mineral Density"),
  list(id = "ieu-b-109", label = "HDL", sheet_name = "HDL Cholesterol"),
  list(id = "ieu-b-110", label = "LDL", sheet_name = "LDL Cholesterol"),
  list(id = "ieu-b-111", label = "Triglycerides", sheet_name = "Triglycerides"),
  list(id = "ebi-a-GCST90014006", label = "HbA1c", sheet_name = "Glycated Hemoglobin")
)
wb <- createWorkbook()
# Define styles for Excel
method_style <- createStyle(textDecoration = "bold", fontSize = 14)
header_style <- createStyle(textDecoration = "bold", fontSize = 14, halign = "center")
make_table <- function(outcome, table_number) {
  id <- outcome$id
  label <- outcome$label
  sheet_name <- paste0("Tab S", table_number, " ", outcome$sheet_name)
  df_MR <- read.csv(paste0(id, "/results/", label, "_MR_GRAPPLE_FDR_p_0.05_summary.csv"))
  df_MVMR_0.05 <- read.csv(paste0(id, "/results/", label, "_MVMR_GRAPPLE_FDR_p_0.05_summary.csv"))
  df_MVMR_0.1 <- read.csv(paste0(id, "/results/", label, "_MVMR_GRAPPLE_FDR_p_0.1_summary.csv"))
  df_MR_strength <- read.csv(paste0(id, "/results/", label, "_final_MR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_UVMR_strength <- read.csv(paste0(id, "/results/", label, "_selection_UVMR_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_onlyBMI_strength <- read.csv(paste0(id, "/results/", label, "_only_BMI_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_MVMR_0.05_strength <- read.csv(paste0(id, "/results/", label, "_final_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_MVMR_0.1_strength <- read.csv(paste0(id, "/results/", label, "_final_MVMR_GRAPPLE_FDR_p_0.1_MVMR_strength.csv"))
  df_MVMR_0.05_all_strength <- read.csv(paste0(id, "/results/", label, "_unique_traits_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_MVMR_0.05_marginal_strength <- read.csv(paste0(id, "/results/", label, "_selection_marginal_p_0.05_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  df_MVMR_0.05_literature_strength <- read.csv(paste0(id, "/results/", label, "_selection_literature_MVMR_GRAPPLE_FDR_p_0.05_MVMR_strength.csv"))
  # Order:
  # 1. UVMR
  # 2. only BMI
  # 3. MVMR-Stepwise F-stats filtering with FDR 0.05
  # 4. MVMR-Stepwise F-stats filtering with FDR 0.1
  # 5. MVMR-Stepwise F-stats filtering, no BMI adjusted, with FDR 0.05
  # 6. MVMR-no F-stats filtering with FDR 0.05
  # 7. MVMR-marginal
  # 8. MVMR-literature
  method_order <- df_MVMR_0.05$method %>% unique()
  method_order <- c("GRAPPLE_5e-08", setdiff(method_order, "GRAPPLE_5e-08"))
  res <- df_MVMR_0.05 %>%
    filter(type == "UVMR") %>%
    select(exposure, b, se, pvalue, method) %>%
    left_join(df_UVMR_strength[, c("F.statistic", "id")], by = c("exposure" = "id")) %>%
    full_join(
      df_MVMR_0.05 %>%
        filter(type == "only_BMI") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_onlyBMI_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MVMR_0.05 %>%
        filter(type == "Stepwise") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MVMR_0.05_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MVMR_0.1 %>%
        filter(type == "Stepwise") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MVMR_0.1_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MR %>%
        filter(type == "Stepwise") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MR_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MVMR_0.05 %>%
        filter(type == "unique_traits") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MVMR_0.05_all_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MVMR_0.05 %>%
        filter(type == "marginal_p_0.05") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MVMR_0.05_marginal_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    full_join(
      df_MVMR_0.05 %>%
        filter(type == "literature") %>%
        select(exposure, b, se, pvalue, method) %>%
        left_join(df_MVMR_0.05_literature_strength[, c("F.statistic", "id")], by = c("exposure" = "id")),
      by = c("exposure", "method")
    ) %>%
    rename("ID" = "exposure") %>%
    mutate(trait = gwasinfo(ID)$trait) %>%
    filter(method != "IVW_T_5e-08" & method != "ESMR_5e-08" & method != "ESMR_optimize_5e-08") %>%
    select(method, ID, trait, everything()) %>%
    mutate(method = factor(method, levels = method_order)) %>%
    arrange(method, desc(ID == "ebi-a-GCST90029070"))
  res <- res %>% droplevels()
  res_list <- list()
  # For each unique method, add a separator row and the corresponding data
  for(m in levels(res$method)) {
    single_method <- res %>% filter(method == m) %>% select(-method)
    method_row <- data.frame(matrix(NA, ncol = ncol(single_method), nrow = 1))
    colnames(method_row) <- colnames(single_method)
    method_row$ID <- m
    combined <- bind_rows(method_row, single_method)
    res_list[[m]] <- combined
  }
  res_final <- bind_rows(res_list)
  colnames(res_final) <- c("ID", "trait", rep(c("beta", "se", "pvalue", "F_stats"),8))
  extra_header <- c("ID", "trait",
                    rep(c("UVMR"),4),
                    rep(c("Only BMI"),4),
                    rep(c("MVMR - main results"),4),
                    rep(c("MVMR - FDR 0.1"),4),
                    rep(c("MVMR - no BMI by default"),4),
                    rep(c("MVMR - MVMR - no F-stats filtering"),4),
                    rep(c("MVMR - Marginal Selection"),4),
                    rep(c("MVMR - Literature"),4))
  # Add worksheet for the current outcome with the specified sheet name
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = as.data.frame(t(extra_header)), startRow = 1, colNames = FALSE)
  writeData(wb, sheet = sheet_name, x = res_final, startRow = 2, colNames = TRUE)
  addStyle(wb, sheet = sheet_name, style = header_style, rows = 1, cols = 3:34, gridExpand = TRUE)
  addStyle(wb, sheet = sheet_name, style = header_style, rows = 1, cols = 1:2, gridExpand = TRUE)

  mergeCells(wb, sheet_name, cols = 3:6, rows = 1)
  mergeCells(wb, sheet_name, cols = 7:10, rows = 1)
  mergeCells(wb, sheet_name, cols = 11:14, rows = 1)
  mergeCells(wb, sheet_name, cols = 15:18, rows = 1)
  mergeCells(wb, sheet_name, cols = 19:22, rows = 1)
  mergeCells(wb, sheet_name, cols = 23:26, rows = 1)
  mergeCells(wb, sheet_name, cols = 27:30, rows = 1)
  mergeCells(wb, sheet_name, cols = 31:34, rows = 1)

  method_rows <- which(!is.na(res_final$ID) & res_final$ID %in% levels(res$method))
  for (row in method_rows) {
    addStyle(wb, sheet = sheet_name, style = method_style, rows = row + 2, cols = 1:ncol(res_final), gridExpand = TRUE)
  }
}
table_number <- 2
for (outcome in outcomes) {
  make_table(outcome, table_number)
  table_number <- table_number + 1
}
# Add Z diff test
res_diff <- read.csv("res_diff.csv")
addWorksheet(wb, "Tab S18 Diff-test")
writeData(wb, sheet = "Tab S18 Diff-test", x = res_diff, startRow = 1, colNames = TRUE)
# Define header style
header_style_S18 <- createStyle(textDecoration = "bold", fontSize = 12, halign = "center")
addStyle(wb, sheet = "Tab S18 Diff-test", style = header_style_S18, rows = 1, cols = 1:ncol(res_diff), gridExpand = TRUE)
# Define bold style for significant p-values
bold_sig_style <- createStyle(textDecoration = "bold", fontColour = "#D7263D")
sig_rows <- which(res_diff$p_diff < 0.05) + 1
if (length(sig_rows) > 0) {
  addStyle(wb, sheet = "Tab S18 Diff-test", style = bold_sig_style, rows = sig_rows, cols = ncol(res_diff), gridExpand = FALSE)
}
# Save the workbook with all outcome results
saveWorkbook(wb, "all_supp_new.xlsx", overwrite = TRUE)


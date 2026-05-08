# Purpose: Baseline diabetes concordance for complete cases.
# Steps: Build 2x2 table, run Fisher test and logistic model, format cell outputs.
# Output: Writes psdcan02_dm concordance table cells.csv.
rm(list=ls());gc();source(".Rprofile")

library(dplyr)

dyads <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS"))

# Baseline diabetes concordance using dyads
baseline_couples <- dyads %>%
  dplyr::select(hhid, dm_biomarker0_wife, dm_biomarker0_husb) %>%
  dplyr::filter(!is.na(dm_biomarker0_wife), !is.na(dm_biomarker0_husb))

# -----------------------------------------------------------------------------
# 1) Frequency table: wife * husband baseline diabetes
# -----------------------------------------------------------------------------
dm_xtab <- with(baseline_couples, table(dm_biomarker0_wife, dm_biomarker0_husb))
print(dm_xtab)

# Relative risk / odds ratio analogue for 2x2 table
# fisher.test returns an exact odds ratio and confidence interval.
dm_fisher <- fisher.test(dm_xtab)
print(dm_fisher)

# -----------------------------------------------------------------------------
# 2) Logistic model: husband baseline DM ~ wife baseline DM
# -----------------------------------------------------------------------------
dm_logit <- glm(dm_biomarker0_husb ~ dm_biomarker0_wife, data = baseline_couples, family = binomial())

coef_est <- coef(dm_logit)["dm_biomarker0_wife"]
coef_se <- sqrt(vcov(dm_logit)["dm_biomarker0_wife", "dm_biomarker0_wife"])
coef_p <- summary(dm_logit)$coefficients["dm_biomarker0_wife", "Pr(>|z|)"]

logit_or_tbl <- tibble::tibble(
  term = "dm_biomarker0_wife (1 vs 0)",
  estimate_log_odds = coef_est,
  odds_ratio = exp(coef_est),
  conf_low_95 = exp(coef_est - 1.96 * coef_se),
  conf_high_95 = exp(coef_est + 1.96 * coef_se),
  p_value = coef_p
)

print(summary(dm_logit))
print(logit_or_tbl)

# -----------------------------------------------------------------------------
# 3) Table for baseline diabetes concordance
# -----------------------------------------------------------------------------

dm_table_labeled <- with(
  baseline_couples,
  table(
    women = factor(dm_biomarker0_wife, levels = c(0, 1), labels = c("No diabetes", "Diabetes")),
    men = factor(dm_biomarker0_husb, levels = c(0, 1), labels = c("No diabetes", "Diabetes"))
  )
)

dm_row_pct <- 100 * prop.table(dm_table_labeled, margin = 1)
dm_col_pct <- 100 * prop.table(dm_table_labeled, margin = 2)

dm_cell_df <- as.data.frame(as.table(dm_table_labeled), stringsAsFactors = FALSE) %>%
  dplyr::rename(Women_Status = women, Men_Status = men, Count = Freq) %>%
  dplyr::mutate(
    overall_pct = 100 * Count / sum(Count),
    row_pct = dm_row_pct[cbind(as.character(Women_Status), as.character(Men_Status))],
    col_pct = dm_col_pct[cbind(as.character(Women_Status), as.character(Men_Status))],
    cell_text = sprintf(
      "%d (%.1f%%)\n(%.1f%%)\n(%.1f%%)",
      Count,
      overall_pct,
      row_pct,
      col_pct
    )
  )

dm_table_display <- dm_cell_df %>%
  dplyr::select(Women_Status, Men_Status, cell_text) %>%
  tidyr::pivot_wider(names_from = Men_Status, values_from = cell_text)

print(dm_table_display)
utils::write.csv(
  dm_cell_df,
  file = "cca/analysis/psdcan02_dm concordance table cells.csv",
  row.names = FALSE
)


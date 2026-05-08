# Purpose: Table 3 complete-case follow-up risk by partner baseline status.
# Steps: Build sex-stratified risk tables and format percentages.
# Output: Writes psdcan03_followup_risk_by_spouse_baseline.csv.
rm(list=ls());gc();source(".Rprofile")

library(dplyr)
library(tidyr)

dyads <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS"))

# ============================================================================
# Table 3: Risk of incident diabetes over follow-up by partner baseline status
# among individuals without diabetes at baseline.
# ============================================================================

build_table_section <- function(data, baseline_label) {
  data %>%
    dplyr::mutate(
      partner_status = factor(
        partner_status,
        levels = c("No diabetes", "Diabetes")
      ),
      outcome_status = factor(
        outcome_status,
        levels = c("No diabetes", "Diabetes")
      )
    ) %>%
    dplyr::count(partner_status, outcome_status, name = "n") %>%
    dplyr::group_by(partner_status) %>%
    dplyr::mutate(
      partner_n = sum(n),
      pct = 100 * n / partner_n
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      baseline_group = baseline_label,
      partner_group = paste0(as.character(partner_status), " (n=", partner_n, ")")
    ) %>%
    dplyr::select(baseline_group, partner_group, outcome_status, pct) %>%
    tidyr::pivot_wider(
      names_from = outcome_status,
      values_from = pct
    ) %>%
    dplyr::rename(
      no_diabetes_pct = `No diabetes`,
      diabetes_pct = Diabetes
    ) %>%
    dplyr::mutate(
      no_diabetes_pct = sprintf("%.1f", no_diabetes_pct),
      diabetes_pct = sprintf("%.1f", diabetes_pct)
    )
}

# Women without diabetes at baseline: follow-up diabetes by husband's baseline status.
women_df <- dyads %>%
  dplyr::filter(dm_biomarker0_wife == 0) %>%
  dplyr::mutate(
    partner_status = dplyr::case_when(
      dm_biomarker0_husb == 0 ~ "No diabetes",
      dm_biomarker0_husb == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    ),
    outcome_status = dplyr::case_when(
      event_DMbiomarker_wife == 0 ~ "No diabetes",
      event_DMbiomarker_wife == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

women_label <- paste0(
  "Women without diabetes at baseline (n=",
  nrow(women_df),
  ")"
)

women_section <- build_table_section(women_df, women_label)

# Men without diabetes at baseline: follow-up diabetes by wife's baseline status.
men_df <- dyads %>%
  dplyr::filter(dm_biomarker0_husb == 0) %>%
  dplyr::mutate(
    partner_status = dplyr::case_when(
      dm_biomarker0_wife == 0 ~ "No diabetes",
      dm_biomarker0_wife == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    ),
    outcome_status = dplyr::case_when(
      event_DMbiomarker_husb == 0 ~ "No diabetes",
      event_DMbiomarker_husb == 1 ~ "Diabetes",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

men_label <- paste0(
  "Men without diabetes at baseline (n=",
  nrow(men_df),
  ")"
)

men_section <- build_table_section(men_df, men_label)

table3_df <- dplyr::bind_rows(women_section, men_section) %>%
  dplyr::rename(
    `Individual's status at baseline` = baseline_group,
    `Individuals with following partner's status at baseline` = partner_group,
    `No diabetes, %` = no_diabetes_pct,
    `Diabetes, %` = diabetes_pct
  )

print(table3_df)

utils::write.csv(
  table3_df,
  file = "cca/analysis/psdcan03_followup_risk_by_spouse_baseline.csv",
  row.names = FALSE
)

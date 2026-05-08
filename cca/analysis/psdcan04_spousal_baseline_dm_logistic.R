# Purpose: Complete-case logistic models for Table 4 (incident DM vs partner status).
# Steps: Fit sex-specific models and format OR tables and plot data.
# Output: Writes Table 4 outputs to cca/analysis/.
rm(list = ls()); gc(); source(".Rprofile")

library(dplyr)
library(ggplot2)

# ============================================================================
# Table 4: Risk of incident diabetes by spousal baseline diabetes status
# among individuals without diabetes at baseline.
# ============================================================================

dyads <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS")
)

extract_partner_or <- function(model_obj, term_name = "partner_dm_statusDiabetes") {
  coef_est <- stats::coef(model_obj)[term_name]
  coef_se <- sqrt(stats::vcov(model_obj)[term_name, term_name])

  tibble::tibble(
    or = exp(coef_est),
    low = exp(coef_est - 1.96 * coef_se),
    high = exp(coef_est + 1.96 * coef_se)
  )
}

fmt_or <- function(or, low, high) {
  sprintf("%.2f (%.2f - %.2f)", or, low, high)
}

fit_sex_logit_models <- function(df, sex_label, include_behavior_terms = FALSE) {
  df_model <- df %>%
    dplyr::filter(
      baseline_dm == 0,
      !is.na(partner_baseline_dm),
      !is.na(dm_incident)
    ) %>%
    dplyr::mutate(
      dm_incident = as.integer(dm_incident),
      partner_dm_status = factor(
        partner_baseline_dm,
        levels = c(0, 1),
        labels = c("No diabetes", "Diabetes")
      ),
      edu_category = as.factor(edu_category),
      famhx_dm = as.factor(famhx_dm),
      morbidity_category = as.factor(morbidity_category),
      bmi_category = as.factor(bmi_category),
      alc_overall = as.factor(alc_overall),
      smk_overall = as.factor(smk_overall)
    )

  model_age <- stats::glm(
    dm_incident ~ partner_dm_status + age,
    data = df_model,
    family = stats::binomial()
  )

  model_1 <- stats::glm(
    dm_incident ~ partner_dm_status + age + edu_category + famhx_dm,
    data = df_model,
    family = stats::binomial()
  )

  if (include_behavior_terms) {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category + alc_overall + smk_overall,
      data = df_model,
      family = stats::binomial()
    )
  } else {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category,
      data = df_model,
      family = stats::binomial()
    )
  }

  age_est <- extract_partner_or(model_age)
  m1_est <- extract_partner_or(model_1)
  m2_est <- extract_partner_or(model_2)

  partner_n <- df_model %>%
    dplyr::count(partner_dm_status, name = "n")

  no_dm_n <- partner_n$n[partner_n$partner_dm_status == "No diabetes"]
  dm_n <- partner_n$n[partner_n$partner_dm_status == "Diabetes"]

  baseline_label <- paste0(
    sex_label,
    " without diabetes at baseline (n=",
    nrow(df_model),
    ")"
  )

  table_rows <- tibble::tibble(
    `Individual's status at baseline` = c(baseline_label, baseline_label),
    `Individuals with following partner's status at baseline` = c(
      paste0("No diabetes (n=", no_dm_n, ")"),
      paste0("Diabetes (n=", dm_n, ")")
    ),
    `Age adjusted OR (95% CI)` = c("1.00", fmt_or(age_est$or, age_est$low, age_est$high)),
    `Model 1 adjusted OR (95% CI)` = c("1.00", fmt_or(m1_est$or, m1_est$low, m1_est$high)),
    `Model 2 adjusted OR (95% CI)` = c("1.00", fmt_or(m2_est$or, m2_est$low, m2_est$high))
  )

  plot_rows <- tibble::tibble(
    sex = sex_label,
    model = c("Age adjusted", "Model 1", "Model 2"),
    or = c(age_est$or, m1_est$or, m2_est$or),
    low = c(age_est$low, m1_est$low, m2_est$low),
    high = c(age_est$high, m1_est$high, m2_est$high)
  )

  list(table_rows = table_rows, plot_rows = plot_rows)
}

women_df <- dyads %>%
  dplyr::transmute(
    baseline_dm = dm_biomarker0_wife,
    partner_baseline_dm = dm_biomarker0_husb,
    dm_incident = event_DMbiomarker_wife,
    age = age_wife,
    edu_category = educcat_wife,
    famhx_dm = famhx_dm_wife,
    morbidity_category = multimorbiditycat_wife,
    bmi_category = bmicat_wife,
    alc_overall = NA,
    smk_overall = NA
  )

men_df <- dyads %>%
  dplyr::transmute(
    baseline_dm = dm_biomarker0_husb,
    partner_baseline_dm = dm_biomarker0_wife,
    dm_incident = event_DMbiomarker_husb,
    age = age_husb,
    edu_category = educcat_husb,
    famhx_dm = famhx_dm_husb,
    morbidity_category = multimorbiditycat_husb,
    bmi_category = bmicat_husb,
    alc_overall = alc_overall_husb,
    smk_overall = smk_overall_husb
  )

women_results <- fit_sex_logit_models(
  df = women_df,
  sex_label = "Women",
  include_behavior_terms = FALSE
)

men_results <- fit_sex_logit_models(
  df = men_df,
  sex_label = "Men",
  include_behavior_terms = TRUE
)

table4_df <- dplyr::bind_rows(
  women_results$table_rows,
  men_results$table_rows
)

utils::write.csv(
  table4_df,
  file = "cca/analysis/psdcan04_table4_spousal_baseline_dm_logistic.csv",
  row.names = FALSE
)

# Purpose: Complete-case Poisson RR models for Table 5.
# Steps: Fit sex-specific models, compute RR with robust SE, format tables.
# Output: Writes Table 5 outputs to cca/analysis/.
rm(list = ls()); gc(); source(".Rprofile")

library(ggplot2)

# ============================================================================
# Table 5: Risk ratio of incident diabetes by spousal baseline diabetes status
# among individuals without diabetes at baseline.
# ============================================================================

dyads <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS")
)

get_model_vcov <- function(model_obj) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    return(sandwich::vcovHC(model_obj, type = "HC0"))
  }

  warning("Package 'sandwich' is not installed. Using model-based variance instead of robust variance.")
  stats::vcov(model_obj)
}

extract_partner_rr <- function(model_obj, term_name = "partner_dm_statusDiabetes") {
  coef_names <- names(stats::coef(model_obj))
  if (!term_name %in% coef_names) {
    stop(paste0("Term not found in model coefficients: ", term_name))
  }

  vcov_mat <- get_model_vcov(model_obj)
  coef_est <- stats::coef(model_obj)[term_name]
  coef_se <- sqrt(vcov_mat[term_name, term_name])

  tibble::tibble(
    rr = exp(coef_est),
    low = exp(coef_est - 1.96 * coef_se),
    high = exp(coef_est + 1.96 * coef_se)
  )
}

fmt_rr <- function(rr, low, high) {
  sprintf("%.2f (%.2f - %.2f)", rr, low, high)
}

fit_sex_rr_models <- function(df, sex_label, include_behavior_terms = FALSE) {
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
    family = stats::poisson(link = "log")
  )

  model_1 <- stats::glm(
    dm_incident ~ partner_dm_status + age + edu_category + famhx_dm,
    data = df_model,
    family = stats::poisson(link = "log")
  )

  if (include_behavior_terms) {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category + alc_overall + smk_overall,
      data = df_model,
      family = stats::poisson(link = "log")
    )
  } else {
    model_2 <- stats::glm(
      dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
        morbidity_category + bmi_category,
      data = df_model,
      family = stats::poisson(link = "log")
    )
  }

  age_est <- extract_partner_rr(model_age)
  m1_est <- extract_partner_rr(model_1)
  m2_est <- extract_partner_rr(model_2)

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
    `Age adjusted RR (95% CI)` = c("1.00", fmt_rr(age_est$rr, age_est$low, age_est$high)),
    `Model 1 adjusted RR (95% CI)` = c("1.00", fmt_rr(m1_est$rr, m1_est$low, m1_est$high)),
    `Model 2 adjusted RR (95% CI)` = c("1.00", fmt_rr(m2_est$rr, m2_est$low, m2_est$high))
  )

  plot_rows <- tibble::tibble(
    sex = sex_label,
    model = c("Age adjusted", "Model 1", "Model 2"),
    rr = c(age_est$rr, m1_est$rr, m2_est$rr),
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

women_results <- fit_sex_rr_models(
  df = women_df,
  sex_label = "Women",
  include_behavior_terms = FALSE
)

men_results <- fit_sex_rr_models(
  df = men_df,
  sex_label = "Men",
  include_behavior_terms = TRUE
)

table5_df <- dplyr::bind_rows(
  women_results$table_rows,
  men_results$table_rows
)


utils::write.csv(
  table5_df,
  file = "cca/analysis/psdcan05_table5_spousal_baseline_dm_rr.csv",
  row.names = FALSE
)

# Archive: Table 6 Cox models for spousal baseline DM using wide complete-case data.
# Steps: Build baseline dyads, fit sex-specific Cox models, format HR rows.
# Output: Prepares Table 6 HR outputs for export.
rm(list=ls());gc();source(".Rprofile")

analytic_df_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre04_analytic dataset wide.RDS"))

analytic_df_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre05_analytic dataset wide.RDS"))

# ============================================================================
# Table 6: Risk of diabetes diagnosis by spousal baseline diabetes status
# among individuals without diabetes at baseline.
# Model 1 adjusts for age, education, family history of diabetes.
# Model 2 adds morbidity category and BMI category.
# ============================================================================

# In this code, morbidity_category is used as multimorbidity_category.

baseline_dyads <- analytic_df_wide %>%
  dplyr::filter(fup == 0) %>%
  dplyr::select(
    hhid,
    carrs,
    female_pid,
    male_pid,
    female_baseline_dm,
    male_baseline_dm,
    female_age,
    male_age,
    female_edu_category,
    male_edu_category,
    female_famhx_dm,
    male_famhx_dm,
    female_morbidity_category,
    male_morbidity_category,
    female_bmi_category,
    male_bmi_category
  ) %>%
  dplyr::distinct()

followup_outcome <- dm_event_or_censor_df %>%
  dplyr::select(carrs, pid, dm_incident, time_to_dm) %>%
  dplyr::distinct(carrs, pid, .keep_all = TRUE)

extract_partner_hr <- function(model_obj, term_name = "partner_dm_statusDiabetes") {
  sm <- summary(model_obj)
  tibble::tibble(
    hr = sm$conf.int[term_name, "exp(coef)"],
    low = sm$conf.int[term_name, "lower .95"],
    high = sm$conf.int[term_name, "upper .95"]
  )
}

fmt_hr <- function(hr, low, high) {
  sprintf("%.2f (%.2f - %.2f)", hr, low, high)
}

fit_sex_models <- function(df, baseline_group_label) {
  df_model <- df %>%
    dplyr::filter(
      baseline_dm == 0,
      !is.na(partner_baseline_dm),
      !is.na(dm_incident),
      !is.na(time_to_dm),
      time_to_dm >= 0
    ) %>%
    dplyr::mutate(
      partner_dm_status = factor(
        partner_baseline_dm,
        levels = c(0, 1),
        labels = c("No diabetes", "Diabetes")
      ),
      edu_category = as.factor(edu_category),
      famhx_dm = as.factor(famhx_dm),
      morbidity_category = as.factor(morbidity_category),
      bmi_category = as.factor(bmi_category)
    )

  # Age-adjusted model
  model_age <- survival::coxph(
    survival::Surv(time_to_dm, dm_incident) ~ partner_dm_status + age,
    data = df_model,
    ties = "efron"
  )

  # Model 1: age + education + family history
  model_1 <- survival::coxph(
    survival::Surv(time_to_dm, dm_incident) ~
      partner_dm_status + age + edu_category + famhx_dm,
    data = df_model,
    ties = "efron"
  )

  # Model 2: Model 1 + morbidity category + BMI category
  model_2 <- survival::coxph(
    survival::Surv(time_to_dm, dm_incident) ~
      partner_dm_status + age + edu_category + famhx_dm +
      morbidity_category + bmi_category,
    data = df_model,
    ties = "efron"
  )

  age_est <- extract_partner_hr(model_age)
  m1_est <- extract_partner_hr(model_1)
  m2_est <- extract_partner_hr(model_2)

  partner_n <- df_model %>%
    dplyr::count(partner_dm_status, name = "n")

  no_dm_n <- partner_n$n[partner_n$partner_dm_status == "No diabetes"]
  dm_n <- partner_n$n[partner_n$partner_dm_status == "Diabetes"]

  tibble::tibble(
    `Individual's status at baseline` = c(
      baseline_group_label,
      baseline_group_label
    ),
    `Individuals with following partner's status at baseline` = c(
      paste0("No diabetes (n=", no_dm_n, ")"),
      paste0("Diabetes (n=", dm_n, ")")
    ),
    `Age adjusted HR (95% CI)` = c(
      "1.00",
      fmt_hr(age_est$hr, age_est$low, age_est$high)
    ),
    `Model 1 adjusted HR (95% CI)` = c(
      "1.00",
      fmt_hr(m1_est$hr, m1_est$low, m1_est$high)
    ),
    `Model 2 adjusted HR (95% CI)` = c(
      "1.00",
      fmt_hr(m2_est$hr, m2_est$low, m2_est$high)
    )
  )
}

women_df <- baseline_dyads %>%
  dplyr::transmute(
    hhid,
    carrs,
    pid = female_pid,
    baseline_dm = female_baseline_dm,
    partner_baseline_dm = male_baseline_dm,
    age = female_age,
    edu_category = female_edu_category,
    famhx_dm = female_famhx_dm,
    morbidity_category = female_morbidity_category,
    bmi_category = female_bmi_category
  ) %>%
  dplyr::left_join(followup_outcome, by = c("carrs", "pid"))

men_df <- baseline_dyads %>%
  dplyr::transmute(
    hhid,
    carrs,
    pid = male_pid,
    baseline_dm = male_baseline_dm,
    partner_baseline_dm = female_baseline_dm,
    age = male_age,
    edu_category = male_edu_category,
    famhx_dm = male_famhx_dm,
    morbidity_category = male_morbidity_category,
    bmi_category = male_bmi_category
  ) %>%
  dplyr::left_join(followup_outcome, by = c("carrs", "pid"))

women_label <- paste0(
  "Women without diabetes at baseline (n=",
  sum(women_df$baseline_dm == 0, na.rm = TRUE),
  ")"
)

men_label <- paste0(
  "Men without diabetes at baseline (n=",
  sum(men_df$baseline_dm == 0, na.rm = TRUE),
  ")"
)

table6_women <- fit_sex_models(women_df, women_label)
table6_men <- fit_sex_models(men_df, men_label)

table6_df <- dplyr::bind_rows(table6_women, table6_men)

print(table6_df)

utils::write.csv(
  table6_df,
  file = "cca/analysis/psdcan05_table6_spousal_baseline_dm_cox.csv",
  row.names = FALSE
)

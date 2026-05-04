rm(list = ls()); gc(); source(".Rprofile")

library(dplyr)

# Debug a single imputation for Table 4 logistic models

dyads_list <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS")
)

dyads_list <- dyads_list[!vapply(dyads_list, is.null, logical(1))]

imp_num <- 1
if (imp_num > length(dyads_list)) stop("imp_num exceeds available imputations")

dyads <- dyads_list[[imp_num]]

fit_one <- function(sex_label, include_behavior_terms = FALSE) {
  if (sex_label == "Women") {
    df <- dyads %>%
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
  } else {
    df <- dyads %>%
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
  }

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

  cat("\n===", sex_label, "===\n")
  print(df_model %>% dplyr::count(partner_dm_status, name = "n"))

  f_age <- dm_incident ~ partner_dm_status + age
  f_m1 <- dm_incident ~ partner_dm_status + age + edu_category + famhx_dm
  f_m2 <- if (include_behavior_terms) {
    dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
      morbidity_category + bmi_category + alc_overall + smk_overall
  } else {
    dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
      morbidity_category + bmi_category
  }

  for (f in list(f_age, f_m1, f_m2)) {
    mf <- stats::model.frame(f, data = df_model, na.action = stats::na.omit)
    cat("\nFormula:", deparse(f), "\n")
    if (!("partner_dm_status" %in% names(mf))) {
      cat("partner_dm_status missing after NA omit\n")
      next
    }
    print(table(mf$partner_dm_status))
    if (length(unique(mf$partner_dm_status)) < 2) {
      cat("No variation in partner_dm_status\n")
      next
    }
    fit <- stats::glm(f, data = mf, family = stats::binomial())
    print(summary(fit)$coefficients)
  }
}

fit_one("Women", include_behavior_terms = FALSE)
fit_one("Men", include_behavior_terms = TRUE)

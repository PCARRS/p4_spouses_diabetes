# Purpose: Complete-case Cox models for spouse status and incident diabetes.
# Steps: Fit sex-specific Cox models and export spouse HRs.
# Output: Writes psdcan06_spouse_cox_models_spouse_hr.csv.
rm(list=ls());gc();source(".Rprofile")

library(dplyr)
library(survival)
library(tibble)

# Load dyad dataset

dyads <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS"))

# Helper to fit Cox model and return spouse-status HR
fit_cox <- function(data, time_var, event_var, spouse_var, covars, model_name) {
  df <- data %>%
    dplyr::filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))

  fml <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ",
           paste(c(spouse_var, covars), collapse = " + "))
  )

  fit <- coxph(fml, data = df)
  summ <- summary(fit)

  if (!(spouse_var %in% rownames(summ$coefficients))) {
    hr <- NA_real_
    lcl <- NA_real_
    ucl <- NA_real_
    pval <- NA_real_
  } else {
    coef_row <- summ$coefficients[spouse_var, , drop = FALSE]
    ci_row <- summ$conf.int[spouse_var, , drop = FALSE]
    hr <- as.numeric(ci_row[, "exp(coef)"])
    lcl <- as.numeric(ci_row[, "lower .95"])
    ucl <- as.numeric(ci_row[, "upper .95"])
    pval <- as.numeric(coef_row[, "Pr(>|z|)"])
  }

  tibble(
    model = model_name,
    n = nrow(df),
    events = sum(df[[event_var]] == 1, na.rm = TRUE),
    hr_spouse = hr,
    hr_lcl = lcl,
    hr_ucl = ucl,
    p_value = pval
  )
}

# Female models
female_base <- dyads %>%
  dplyr::filter(dm_biomarker0_wife == 0, dyad_at_risk == 1)

female_models <- bind_rows(
  fit_cox(
    female_base,
    time_var = "time_to_dm_wife",
    event_var = "event_DMbiomarker_wife",
    spouse_var = "dm_biomarker0_husb",
    covars = c(),
    model_name = "Female: spouse status"
  ),
  fit_cox(
    female_base,
    time_var = "time_to_dm_wife",
    event_var = "event_DMbiomarker_wife",
    spouse_var = "dm_biomarker0_husb",
    covars = c("age_wife"),
    model_name = "Female: spouse status + age"
  ),
  fit_cox(
    female_base,
    time_var = "time_to_dm_wife",
    event_var = "event_DMbiomarker_wife",
    spouse_var = "dm_biomarker0_husb",
    covars = c("age_wife", "educcat_wife", "famhx_dm_wife"),
    model_name = "Female: spouse status + Model 1"
  ),
  fit_cox(
    female_base,
    time_var = "time_to_dm_wife",
    event_var = "event_DMbiomarker_wife",
    spouse_var = "dm_biomarker0_husb",
    covars = c("age_wife", "educcat_wife", "famhx_dm_wife",
               "multimorbiditycat_wife", "bmicat_wife"),
    model_name = "Female: spouse status + Model 2"
  )
)

# Male models
male_base <- dyads %>%
  dplyr::filter(dm_biomarker0_husb == 0, dyad_at_risk == 1)

male_models <- bind_rows(
  fit_cox(
    male_base,
    time_var = "time_to_dm_husb",
    event_var = "event_DMbiomarker_husb",
    spouse_var = "dm_biomarker0_wife",
    covars = c(),
    model_name = "Male: spouse status"
  ),
  fit_cox(
    male_base,
    time_var = "time_to_dm_husb",
    event_var = "event_DMbiomarker_husb",
    spouse_var = "dm_biomarker0_wife",
    covars = c("age_husb"),
    model_name = "Male: spouse status + age"
  ),
  fit_cox(
    male_base,
    time_var = "time_to_dm_husb",
    event_var = "event_DMbiomarker_husb",
    spouse_var = "dm_biomarker0_wife",
    covars = c("age_husb", "educcat_husb", "famhx_dm_husb"),
    model_name = "Male: spouse status + Model 1"
  ),
  fit_cox(
    male_base,
    time_var = "time_to_dm_husb",
    event_var = "event_DMbiomarker_husb",
    spouse_var = "dm_biomarker0_wife",
    covars = c("age_husb", "educcat_husb", "famhx_dm_husb",
               "smk_overall_husb", "alc_overall_husb",
               "multimorbiditycat_husb", "bmicat_husb"),
    model_name = "Male: spouse status + Model 2"
  )
)

cox_results <- bind_rows(female_models, male_models)

print(cox_results)

write.csv(
  cox_results,
  "cca/analysis/psdcan06_spouse_cox_models_spouse_hr.csv",
  row.names = FALSE
)

# Archive: Cox PH model for incident diabetes using long-format complete-case data.
# Steps: Build event/censor dataset, fit Cox with clustered SE, plot HRs and check PH.
# Output: Saves an HR plot and prints hazard ratio tables.
rm(list=ls());gc();source(".Rprofile")

source("cca/preprocessing/archive/psdcpre03_dm incidence df.R")

# ============================================================================
# Cox PH model for incident diabetes among baseline diabetes-free participants.
# ============================================================================

# One event/censor row per participant from preprocessing step.
event_censor_df <- dm_event_or_censor_df %>%
  dplyr::distinct(carrs, pid, .keep_all = TRUE)

# Baseline covariates at fup = 0 (time-fixed covariates for Cox PH).
baseline_covars <- analytic_df %>%
  dplyr::filter(fup == 0) %>%
  dplyr::select(
    carrs,
    pid,
    hhid,
    sex,
    age,
    bmi,
    edu_category,
    employ_category,
    hhincome,
    site,
    dm_biomarker,
    famhx_dm,
    baseline_dm
  ) %>%
  dplyr::rename(bmi_baseline = bmi) %>%
  dplyr::distinct(carrs, pid, .keep_all = TRUE)

cox_df <- event_censor_df %>%
  dplyr::left_join(baseline_covars, by = c("carrs", "pid", "hhid")) %>%
  dplyr::filter(
    baseline_dm == 0,
    !is.na(time_to_dm),
    !is.na(dm_incident),
    time_to_dm >= 0
  ) %>%
  dplyr::mutate(
    dm_incident = as.integer(dm_incident),
    fup_duration = time_to_dm,
    sex = as.factor(sex),
    edu_category = as.factor(edu_category),
    employ_category = as.factor(employ_category),
    hhincome = as.factor(hhincome),
    site = as.factor(site),
    carrs = as.factor(carrs),
    famhx_dm = as.factor(famhx_dm)
  )

# Main model: robust SE clustered at household level to account for within-couple correlation.
# Note: fup_duration equals the survival time variable (time_to_dm), so it is not included
# on the right-hand side of the Cox model.
cox_model <- survival::coxph(
  survival::Surv(time_to_dm, dm_incident) ~
    sex + age + bmi_baseline + edu_category + employ_category + hhincome + site +
    carrs + famhx_dm + cluster(hhid),
  data = cox_df,
  ties = "efron"
)

cox_summary <- summary(cox_model)
print(cox_summary)

cox_hr_tbl <- tibble::tibble(
  term = rownames(cox_summary$coefficients),
  hazard_ratio = cox_summary$conf.int[, "exp(coef)"],
  conf_low_95 = cox_summary$conf.int[, "lower .95"],
  conf_high_95 = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)
print(cox_hr_tbl)

# Plot hazard ratios with 95% CI.
hr_terms_excluded <- cox_hr_tbl %>%
  dplyr::arrange(hazard_ratio) %>%
  dplyr::slice_head(n = 2) %>%
  dplyr::pull(term)

hr_plot_df <- cox_hr_tbl %>%
  dplyr::filter(!term %in% hr_terms_excluded) %>%
  dplyr::mutate(term = factor(term, levels = rev(term)))

print(hr_terms_excluded)

hr_plot <- ggplot2::ggplot(
  hr_plot_df,
  ggplot2::aes(x = term, y = hazard_ratio, ymin = conf_low_95, ymax = conf_high_95)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  ggplot2::geom_pointrange(size = 0.35) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Cox PH: Hazard Ratios for Incident Diabetes",
    x = "Covariate",
    y = "Hazard Ratio (log scale)"
  ) +
  ggplot2::theme_bw(base_size = 11)

print(hr_plot)

ggplot2::ggsave(
  filename = "cca/analysis/psdcan03_coxph_dm_incidence_hr_plot.png",
  plot = hr_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Proportional hazards assumption check.
cox_ph_test <- survival::cox.zph(cox_model)
print(cox_ph_test)

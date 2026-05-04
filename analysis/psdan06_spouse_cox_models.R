rm(list=ls());gc();source(".Rprofile")

library(survival)

# Cox models pooled across imputations

dyads_list <- readRDS(paste0(path_spouses_diabetes_folder,"/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS"))
dyads_list <- dyads_list[!vapply(dyads_list, is.null, logical(1))]

pool_rubin <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]
  se <- se[ok]
  m <- length(est)
  if (m == 0) return(list(qbar = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
  qbar <- mean(est)
  ubar <- mean(se^2)
  b <- stats::var(est)
  if (is.na(b)) b <- 0
  tvar <- ubar + (1 + 1 / m) * b
  se_p <- sqrt(tvar)
  if (b == 0) {
    lcl <- qbar - 1.96 * se_p
    ucl <- qbar + 1.96 * se_p
  } else {
    df <- (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2
    tcrit <- stats::qt(0.975, df = df)
    lcl <- qbar - tcrit * se_p
    ucl <- qbar + tcrit * se_p
  }
  list(qbar = qbar, se = se_p, lcl = lcl, ucl = ucl)
}

fit_cox_one <- function(dyads, time_var, event_var, spouse_var, covars, model_name) {
  df <- dyads %>%
    dplyr::filter(!is.na(.data[[time_var]]), !is.na(.data[[event_var]]))

  fml <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ",
           paste(c(spouse_var, covars), collapse = " + "))
  )

  fit <- coxph(fml, data = df)
  term <- spouse_var
  est <- stats::coef(fit)[term]
  se <- sqrt(stats::vcov(fit)[term, term])

  list(
    model = model_name,
    est = est,
    se = se,
    n = nrow(df),
    events = sum(df[[event_var]] == 1, na.rm = TRUE)
  )
}

run_models <- function(dyads) {
  female_base <- dyads %>%
    dplyr::filter(dm_biomarker0_wife == 0, dyad_at_risk == 1)

  male_base <- dyads %>%
    dplyr::filter(dm_biomarker0_husb == 0, dyad_at_risk == 1)

  list(
    fit_cox_one(female_base, "time_to_dm_wife", "event_DMbiomarker_wife", "dm_biomarker0_husb", c(), "Female: spouse status"),
    fit_cox_one(female_base, "time_to_dm_wife", "event_DMbiomarker_wife", "dm_biomarker0_husb", c("age_wife"), "Female: spouse status + age"),
    fit_cox_one(female_base, "time_to_dm_wife", "event_DMbiomarker_wife", "dm_biomarker0_husb", c("age_wife", "educcat_wife", "famhx_dm_wife"), "Female: spouse status + Model 1"),
    fit_cox_one(female_base, "time_to_dm_wife", "event_DMbiomarker_wife", "dm_biomarker0_husb", c("age_wife", "educcat_wife", "famhx_dm_wife", "multimorbiditycat_wife", "bmicat_wife"), "Female: spouse status + Model 2"),
    fit_cox_one(male_base, "time_to_dm_husb", "event_DMbiomarker_husb", "dm_biomarker0_wife", c(), "Male: spouse status"),
    fit_cox_one(male_base, "time_to_dm_husb", "event_DMbiomarker_husb", "dm_biomarker0_wife", c("age_husb"), "Male: spouse status + age"),
    fit_cox_one(male_base, "time_to_dm_husb", "event_DMbiomarker_husb", "dm_biomarker0_wife", c("age_husb", "educcat_husb", "famhx_dm_husb"), "Male: spouse status + Model 1"),
    fit_cox_one(male_base, "time_to_dm_husb", "event_DMbiomarker_husb", "dm_biomarker0_wife", c("age_husb", "educcat_husb", "famhx_dm_husb", "smk_overall_husb", "alc_overall_husb", "multimorbiditycat_husb", "bmicat_husb"), "Male: spouse status + Model 2")
  )
}

all_imp <- lapply(dyads_list, run_models)

model_names <- all_imp[[1]]
model_labels <- vapply(model_names, `[[`, character(1), "model")

rows <- lapply(seq_along(model_labels), function(i) {
  ests <- vapply(all_imp, function(x) x[[i]]$est, numeric(1))
  ses <- vapply(all_imp, function(x) x[[i]]$se, numeric(1))
  ns <- vapply(all_imp, function(x) x[[i]]$n, numeric(1))
  evs <- vapply(all_imp, function(x) x[[i]]$events, numeric(1))

  pool <- pool_rubin(ests, ses)
  hr <- exp(pool$qbar)
  lcl <- exp(pool$lcl)
  ucl <- exp(pool$ucl)
  z <- pool$qbar / pool$se
  pval <- 2 * (1 - stats::pnorm(abs(z)))

  tibble::tibble(
    model = model_labels[i],
    n = round(mean(ns, na.rm = TRUE)),
    events = round(mean(evs, na.rm = TRUE)),
    hr_spouse = hr,
    hr_lcl = lcl,
    hr_ucl = ucl,
    p_value = pval
  )
})

cox_results <- dplyr::bind_rows(rows)

write.csv(
  cox_results,
  "analysis/psdan06_spouse_cox_models_spouse_hr.csv",
  row.names = FALSE
)

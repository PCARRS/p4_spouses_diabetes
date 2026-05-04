rm(list = ls()); gc(); source(".Rprofile")

# Table 5 pooled RR models

dyads_list <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS")
)
dyads_list <- dyads_list[!vapply(dyads_list, is.null, logical(1))]

get_model_vcov <- function(model_obj) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    return(sandwich::vcovHC(model_obj, type = "HC0"))
  }
  stats::vcov(model_obj)
}

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

fit_models_one <- function(dyads, sex_label, include_behavior_terms = FALSE, imp_num = NA_integer_) {
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

  fit_log_rr <- function(formula, data) {
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    if (!("partner_dm_status" %in% names(mf))) {
      return(list(est = NA_real_, se = NA_real_, term_present = FALSE))
    }
    if (length(unique(mf$partner_dm_status)) < 2) {
      return(list(est = NA_real_, se = NA_real_, term_present = FALSE))
    }
    fit <- stats::glm(formula, data = mf, family = stats::poisson(link = "log"))
    term <- "partner_dm_statusDiabetes"
    coef_names <- names(stats::coef(fit))
    if (!(term %in% coef_names)) {
      return(list(est = NA_real_, se = NA_real_, term_present = FALSE))
    }
    coef_tbl <- summary(fit)$coefficients
    est <- coef_tbl[term, "Estimate"]
    se <- coef_tbl[term, "Std. Error"]
    if (!is.finite(est)) {
      est <- stats::coef(fit)[term]
    }
    if (!is.finite(se)) {
      vcov_mat <- get_model_vcov(fit)
      se <- sqrt(vcov_mat[term, term])
    }
    list(est = est, se = se, term_present = TRUE)
  }

  partner_n <- df_model %>%
    dplyr::count(partner_dm_status, name = "n")

  age_est <- fit_log_rr(dm_incident ~ partner_dm_status + age, df_model)
  m1_est <- fit_log_rr(dm_incident ~ partner_dm_status + age + edu_category + famhx_dm, df_model)
  m2_est <- if (include_behavior_terms) {
    fit_log_rr(dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
                 morbidity_category + bmi_category + alc_overall + smk_overall, df_model)
  } else {
    fit_log_rr(dm_incident ~ partner_dm_status + age + edu_category + famhx_dm +
                 morbidity_category + bmi_category, df_model)
  }

  n_no_dm <- partner_n$n[partner_n$partner_dm_status == "No diabetes"]
  n_dm <- partner_n$n[partner_n$partner_dm_status == "Diabetes"]

  tibble::tibble(
    imp = imp_num,
    sex = sex_label,
    model = c("Age", "Model1", "Model2"),
    log_rr = c(age_est$est, m1_est$est, m2_est$est),
    log_rr_se = c(age_est$se, m1_est$se, m2_est$se),
    term_present = c(age_est$term_present, m1_est$term_present, m2_est$term_present),
    n_total = nrow(df_model),
    n_no_dm = ifelse(length(n_no_dm) == 0, NA_real_, n_no_dm),
    n_dm = ifelse(length(n_dm) == 0, NA_real_, n_dm)
  )
}

collect_sex <- function(sex_label, include_behavior_terms) {
  res_tbl <- dplyr::bind_rows(lapply(seq_along(dyads_list), function(i) {
    fit_models_one(dyads_list[[i]], sex_label = sex_label, include_behavior_terms = include_behavior_terms, imp_num = i)
  }))

  utils::write.csv(
    res_tbl,
    file = paste0("analysis/psdan05_rr_models_", tolower(sex_label), ".csv"),
    row.names = FALSE
  )

  pool_tbl <- res_tbl %>%
    dplyr::group_by(model) %>%
    dplyr::reframe(
      pool = list(pool_rubin(log_rr, log_rr_se)),
      n_total = round(mean(n_total, na.rm = TRUE)),
      n_no_dm = round(mean(n_no_dm, na.rm = TRUE)),
      n_dm = round(mean(n_dm, na.rm = TRUE))
    ) %>%
    dplyr::mutate(
      log_rr = vapply(pool, `[[`, numeric(1), "qbar"),
      lcl = vapply(pool, `[[`, numeric(1), "lcl"),
      ucl = vapply(pool, `[[`, numeric(1), "ucl")
    )

  pool_tbl
}

women <- collect_sex("Women", include_behavior_terms = FALSE)
men <- collect_sex("Men", include_behavior_terms = TRUE)

fmt_rr <- function(log_rr, lcl, ucl) {
  sprintf("%.2f (%.2f - %.2f)", exp(log_rr), exp(lcl), exp(ucl))
}

build_rows <- function(sex_label, res_tbl) {
  n_total <- unique(res_tbl$n_total)
  n_no_dm <- unique(res_tbl$n_no_dm)
  n_dm <- unique(res_tbl$n_dm)
  baseline_label <- paste0(sex_label, " without diabetes at baseline (n=", n_total, ")")
  age <- res_tbl %>% dplyr::filter(model == "Age")
  m1 <- res_tbl %>% dplyr::filter(model == "Model1")
  m2 <- res_tbl %>% dplyr::filter(model == "Model2")

  tibble::tibble(
    `Individual's status at baseline` = c(baseline_label, baseline_label),
    `Individuals with following partner's status at baseline` = c(
      paste0("No diabetes (n=", n_no_dm, ")"),
      paste0("Diabetes (n=", n_dm, ")")
    ),
    `Age adjusted RR (95% CI)` = c("1.00", fmt_rr(age$log_rr, age$lcl, age$ucl)),
    `Model 1 adjusted RR (95% CI)` = c("1.00", fmt_rr(m1$log_rr, m1$lcl, m1$ucl)),
    `Model 2 adjusted RR (95% CI)` = c("1.00", fmt_rr(m2$log_rr, m2$lcl, m2$ucl))
  )
}

table5_df <- dplyr::bind_rows(
  build_rows("Women", women),
  build_rows("Men", men)
)

utils::write.csv(
  table5_df,
  file = "analysis/psdan05_table5_spousal_baseline_dm_rr.csv",
  row.names = FALSE
)

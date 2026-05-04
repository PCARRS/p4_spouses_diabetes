rm(list=ls());gc();source(".Rprofile")

library(dplyr)
library(tibble)

# Pooled descriptive characteristics across imputed dyads

dyads_list <- readRDS(paste0(path_spouses_diabetes_folder,"/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS"))
dyads_list <- dyads_list[!vapply(dyads_list, is.null, logical(1))]

continuous_vars <- c("age", "bmi", "waist_cm")
proportion_vars <- c("smk_overall", "alc_overall", "famhx_htn", "famhx_dm", "famhx_cvd",
                     "chd", "cva", "ckd", "dm_biomarker0", "htn")
grouped_vars <- c("educcat", "employocccat", "multimorbiditycat", "bmicat")

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

pool_fisher_p <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0]
  if (length(pvals) == 0) return(NA_real_)
  stat <- -2 * sum(log(pvals))
  stats::pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
}

cont_one <- function(dyads) {
  dplyr::bind_rows(lapply(continuous_vars, function(var) {
    female_col <- paste0(var, "_wife")
    male_col <- paste0(var, "_husb")
    if (!(female_col %in% names(dyads)) || !(male_col %in% names(dyads))) {
      return(tibble(variable = var, f_mean = NA_real_, f_sd = NA_real_, f_se = NA_real_,
                    m_mean = NA_real_, m_sd = NA_real_, m_se = NA_real_,
                    r = NA_real_, r_se = NA_real_, f_miss = NA_real_, m_miss = NA_real_))
    }
    f_vals <- dyads[[female_col]]
    m_vals <- dyads[[male_col]]
    f_n <- sum(!is.na(f_vals))
    m_n <- sum(!is.na(m_vals))
    f_mean <- mean(f_vals, na.rm = TRUE)
    m_mean <- mean(m_vals, na.rm = TRUE)
    f_sd <- stats::sd(f_vals, na.rm = TRUE)
    m_sd <- stats::sd(m_vals, na.rm = TRUE)
    f_se <- if (f_n > 0) f_sd / sqrt(f_n) else NA_real_
    m_se <- if (m_n > 0) m_sd / sqrt(m_n) else NA_real_
    r <- suppressWarnings(cor(f_vals, m_vals, use = "pairwise.complete.obs"))
    n_pair <- sum(complete.cases(f_vals, m_vals))
    r_se <- if (is.finite(r) && n_pair > 3) 1 / sqrt(n_pair - 3) else NA_real_
    tibble(
      variable = var,
      f_mean = f_mean,
      f_sd = f_sd,
      f_se = f_se,
      m_mean = m_mean,
      m_sd = m_sd,
      m_se = m_se,
      r = r,
      r_se = r_se,
      f_miss = mean(is.na(f_vals)),
      m_miss = mean(is.na(m_vals))
    )
  }))
}

bin_one <- function(dyads) {
  dplyr::bind_rows(lapply(proportion_vars, function(var) {
    female_col <- paste0(var, "_wife")
    male_col <- paste0(var, "_husb")
    if (!(female_col %in% names(dyads)) || !(male_col %in% names(dyads))) {
      return(tibble(variable = var, f_p = NA_real_, f_se = NA_real_, f_n = NA_real_,
                    m_p = NA_real_, m_se = NA_real_, m_n = NA_real_,
                    log_or = NA_real_, log_or_se = NA_real_, f_miss = NA_real_, m_miss = NA_real_))
    }
    f_vals <- dyads[[female_col]]
    m_vals <- dyads[[male_col]]
    f_n <- sum(!is.na(f_vals))
    m_n <- sum(!is.na(m_vals))
    f_p <- mean(f_vals == 1, na.rm = TRUE)
    m_p <- mean(m_vals == 1, na.rm = TRUE)
    f_se <- if (f_n > 0) sqrt(f_p * (1 - f_p) / f_n) else NA_real_
    m_se <- if (m_n > 0) sqrt(m_p * (1 - m_p) / m_n) else NA_real_
    log_or <- NA_real_
    log_or_se <- NA_real_
    if (f_n > 0 && m_n > 0) {
      glm_fit <- tryCatch(
        stats::glm(m_vals ~ f_vals, family = stats::binomial()),
        error = function(e) NULL
      )
      if (!is.null(glm_fit)) {
        coef_est <- stats::coef(glm_fit)["f_vals"]
        coef_se <- sqrt(stats::vcov(glm_fit)["f_vals", "f_vals"])
        log_or <- as.numeric(coef_est)
        log_or_se <- as.numeric(coef_se)
      }
    }
    tibble(
      variable = var,
      f_p = f_p,
      f_se = f_se,
      f_n = f_n,
      m_p = m_p,
      m_se = m_se,
      m_n = m_n,
      log_or = log_or,
      log_or_se = log_or_se,
      f_miss = mean(is.na(f_vals)),
      m_miss = mean(is.na(m_vals))
    )
  }))
}

cat_one <- function(dyads) {
  dplyr::bind_rows(lapply(grouped_vars, function(var) {
    female_col <- paste0(var, "_wife")
    male_col <- paste0(var, "_husb")
    if (!(female_col %in% names(dyads)) || !(male_col %in% names(dyads))) {
      return(tibble(variable = var, level = NA_character_, f_p = NA_real_, f_se = NA_real_,
                    m_p = NA_real_, m_se = NA_real_, pval = NA_real_, f_miss = NA_real_, m_miss = NA_real_))
    }
    f_vals <- dyads[[female_col]]
    m_vals <- dyads[[male_col]]
    all_levels <- sort(unique(c(f_vals, m_vals)))
    all_levels <- all_levels[!is.na(all_levels)]

    chi_tab <- table(data.frame(
      var = c(f_vals, m_vals),
      sex = rep(c("female", "male"), times = c(length(f_vals), length(m_vals)))
    ))
    chi_result <- tryCatch(stats::chisq.test(chi_tab), error = function(e) NULL)
    pval <- if (!is.null(chi_result)) chi_result$p.value else NA_real_

    f_denom <- sum(!is.na(f_vals))
    m_denom <- sum(!is.na(m_vals))

    dplyr::bind_rows(lapply(seq_along(all_levels), function(i) {
      level <- all_levels[i]
      f_p <- if (f_denom > 0) sum(f_vals == level, na.rm = TRUE) / f_denom else NA_real_
      m_p <- if (m_denom > 0) sum(m_vals == level, na.rm = TRUE) / m_denom else NA_real_
      f_se <- if (f_denom > 0) sqrt(f_p * (1 - f_p) / f_denom) else NA_real_
      m_se <- if (m_denom > 0) sqrt(m_p * (1 - m_p) / m_denom) else NA_real_
      tibble(
        variable = var,
        level = as.character(level),
        f_p = f_p,
        f_se = f_se,
        m_p = m_p,
        m_se = m_se,
        pval = pval,
        f_miss = mean(is.na(f_vals)),
        m_miss = mean(is.na(m_vals))
      )
    }))
  }))
}

cont_all <- dplyr::bind_rows(lapply(dyads_list, cont_one), .id = "imp")
bin_all <- dplyr::bind_rows(lapply(dyads_list, bin_one), .id = "imp")
cat_all <- dplyr::bind_rows(lapply(dyads_list, cat_one), .id = "imp")

continuous_tbl <- cont_all %>%
  dplyr::group_by(variable) %>%
  dplyr::reframe(
    f_pool = list(pool_rubin(f_mean, f_se)),
    m_pool = list(pool_rubin(m_mean, m_se)),
    r_pool = list(pool_rubin(atanh(r), r_se)),
    f_sd = mean(f_sd, na.rm = TRUE),
    m_sd = mean(m_sd, na.rm = TRUE),
    f_miss = mean(f_miss, na.rm = TRUE),
    m_miss = mean(m_miss, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    f_mean = vapply(f_pool, `[[`, numeric(1), "qbar"),
    m_mean = vapply(m_pool, `[[`, numeric(1), "qbar"),
    r = tanh(vapply(r_pool, `[[`, numeric(1), "qbar"))
  ) %>%
  dplyr::transmute(
    variable = variable,
    level = NA_character_,
    female = sprintf("%.1f (%.1f)", f_mean, f_sd),
    male = sprintf("%.1f (%.1f)", m_mean, m_sd),
    compare = sprintf("%.3f", r),
    missing_female = sprintf("%.1f%%", 100 * f_miss),
    missing_male = sprintf("%.1f%%", 100 * m_miss)
  )

binary_tbl <- bin_all %>%
  dplyr::group_by(variable) %>%
  dplyr::reframe(
    f_pool = list(pool_rubin(f_p, f_se)),
    m_pool = list(pool_rubin(m_p, m_se)),
    log_or_pool = list(pool_rubin(log_or, log_or_se)),
    f_n = round(mean(f_n, na.rm = TRUE)),
    m_n = round(mean(m_n, na.rm = TRUE)),
    f_miss = mean(f_miss, na.rm = TRUE),
    m_miss = mean(m_miss, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    f_p = vapply(f_pool, `[[`, numeric(1), "qbar"),
    m_p = vapply(m_pool, `[[`, numeric(1), "qbar"),
    log_or = vapply(log_or_pool, `[[`, numeric(1), "qbar"),
    log_or_se = vapply(log_or_pool, `[[`, numeric(1), "se")
  ) %>%
  dplyr::mutate(
    or = exp(log_or),
    or_lcl = exp(log_or - 1.96 * log_or_se),
    or_ucl = exp(log_or + 1.96 * log_or_se)
  ) %>%
  dplyr::transmute(
    variable = variable,
    level = NA_character_,
    female = sprintf("%d (%.1f%%)", f_n, 100 * f_p),
    male = sprintf("%d (%.1f%%)", m_n, 100 * m_p),
    compare = ifelse(is.finite(or), sprintf("%.2f (%.2f, %.2f)", or, or_lcl, or_ucl), "NA"),
    missing_female = sprintf("%.1f%%", 100 * f_miss),
    missing_male = sprintf("%.1f%%", 100 * m_miss)
  )

categorical_tbl <- cat_all %>%
  dplyr::group_by(variable, level) %>%
  dplyr::reframe(
    f_pool = list(pool_rubin(f_p, f_se)),
    m_pool = list(pool_rubin(m_p, m_se)),
    pval = pool_fisher_p(pval),
    f_miss = mean(f_miss, na.rm = TRUE),
    m_miss = mean(m_miss, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    f_p = vapply(f_pool, `[[`, numeric(1), "qbar"),
    m_p = vapply(m_pool, `[[`, numeric(1), "qbar"),
    pval_fmt = ifelse(is.na(pval), "NA", ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)))
  ) %>%
  dplyr::transmute(
    variable = variable,
    level = level,
    female = sprintf("%.1f%%", 100 * f_p),
    male = sprintf("%.1f%%", 100 * m_p),
    compare = pval_fmt,
    missing_female = sprintf("%.1f%%", 100 * f_miss),
    missing_male = sprintf("%.1f%%", 100 * m_miss)
  ) %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(
    compare = ifelse(row_number() == 1, compare, ""),
    missing_female = ifelse(row_number() == 1, missing_female, ""),
    missing_male = ifelse(row_number() == 1, missing_male, "")
  ) %>%
  dplyr::ungroup()

table1 <- dplyr::bind_rows(continuous_tbl, binary_tbl, categorical_tbl)

write.csv(table1, "analysis/psdan01_descriptive characteristics.csv", row.names = FALSE)

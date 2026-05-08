# Purpose: Estimate baseline diabetes concordance across imputations.
# Steps: Build 2x2 counts, fit logistic models per imputation, pool ORs and p-values.
# Output: Writes concordance table cells and pooled OR outputs to analysis/.
rm(list=ls());gc();source(".Rprofile")

# Baseline diabetes concordance pooled across imputations

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

pool_fisher_p <- function(pvals) {
  pvals <- pvals[is.finite(pvals) & pvals > 0]
  if (length(pvals) == 0) return(NA_real_)
  stat <- -2 * sum(log(pvals))
  stats::pchisq(stat, df = 2 * length(pvals), lower.tail = FALSE)
}

one_imp <- function(dyads) {
  baseline_couples <- dyads %>%
    dplyr::select(hhid, dm_biomarker0_wife, dm_biomarker0_husb) %>%
    dplyr::filter(!is.na(dm_biomarker0_wife), !is.na(dm_biomarker0_husb))

  tab <- with(baseline_couples, table(dm_biomarker0_wife, dm_biomarker0_husb))
  c00 <- tab["0","0"]
  c01 <- tab["0","1"]
  c10 <- tab["1","0"]
  c11 <- tab["1","1"]

  glm_fit <- tryCatch(
    stats::glm(dm_biomarker0_husb ~ dm_biomarker0_wife, data = baseline_couples, family = stats::binomial()),
    error = function(e) NULL
  )

  log_or <- NA_real_
  log_or_se <- NA_real_
  pval <- NA_real_
  if (!is.null(glm_fit)) {
    log_or <- stats::coef(glm_fit)["dm_biomarker0_wife"]
    log_or_se <- sqrt(stats::vcov(glm_fit)["dm_biomarker0_wife", "dm_biomarker0_wife"])
    pval <- summary(glm_fit)$coefficients["dm_biomarker0_wife", "Pr(>|z|)"]
  }

  fisher_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)

  list(
    counts = c(c00 = c00, c01 = c01, c10 = c10, c11 = c11),
    log_or = log_or,
    log_or_se = log_or_se,
    pval = pval,
    fisher_p = fisher_p
  )
}

imp_results <- lapply(dyads_list, one_imp)

counts_mat <- do.call(rbind, lapply(imp_results, function(x) x$counts))
counts_mean <- colMeans(counts_mat, na.rm = TRUE)

log_or_pool <- pool_rubin(
  vapply(imp_results, `[[`, numeric(1), "log_or"),
  vapply(imp_results, `[[`, numeric(1), "log_or_se")
)

log_or <- log_or_pool$qbar
log_or_se <- log_or_pool$se
or <- exp(log_or)
or_lcl <- exp(log_or_pool$lcl)
or_ucl <- exp(log_or_pool$ucl)

p_pooled <- 2 * (1 - stats::pnorm(abs(log_or / log_or_se)))
p_fisher <- pool_fisher_p(vapply(imp_results, `[[`, numeric(1), "fisher_p"))

logit_or_tbl <- tibble::tibble(
  term = "dm_biomarker0_wife (1 vs 0)",
  odds_ratio = or,
  conf_low_95 = or_lcl,
  conf_high_95 = or_ucl,
  p_value = p_pooled,
  p_value_fisher = p_fisher
)

# Build pooled cell table
count_total <- sum(counts_mean)
overall_pct <- 100 * counts_mean / count_total

row_totals <- c(
  row0 = counts_mean["c00"] + counts_mean["c01"],
  row1 = counts_mean["c10"] + counts_mean["c11"]
)
col_totals <- c(
  col0 = counts_mean["c00"] + counts_mean["c10"],
  col1 = counts_mean["c01"] + counts_mean["c11"]
)

row_pct <- c(
  c00 = 100 * counts_mean["c00"] / row_totals["row0"],
  c01 = 100 * counts_mean["c01"] / row_totals["row0"],
  c10 = 100 * counts_mean["c10"] / row_totals["row1"],
  c11 = 100 * counts_mean["c11"] / row_totals["row1"]
)

col_pct <- c(
  c00 = 100 * counts_mean["c00"] / col_totals["col0"],
  c01 = 100 * counts_mean["c01"] / col_totals["col1"],
  c10 = 100 * counts_mean["c10"] / col_totals["col0"],
  c11 = 100 * counts_mean["c11"] / col_totals["col1"]
)

cell_df <- tibble::tibble(
  Women_Status = rep(c("No diabetes", "Diabetes"), each = 2),
  Men_Status = rep(c("No diabetes", "Diabetes"), times = 2),
  Count = c(counts_mean["c00"], counts_mean["c01"], counts_mean["c10"], counts_mean["c11"]),
  overall_pct = c(overall_pct["c00"], overall_pct["c01"], overall_pct["c10"], overall_pct["c11"]),
  row_pct = c(row_pct["c00"], row_pct["c01"], row_pct["c10"], row_pct["c11"]),
  col_pct = c(col_pct["c00"], col_pct["c01"], col_pct["c10"], col_pct["c11"])
) %>%
  dplyr::mutate(
    cell_text = sprintf(
      "%d (%.1f%%)\n(%.1f%%)\n(%.1f%%)",
      round(Count), overall_pct, row_pct, col_pct
    )
  )

utils::write.csv(
  cell_df,
  file = "analysis/psdan02_dm concordance table cells.csv",
  row.names = FALSE
)

utils::write.csv(
  logit_or_tbl,
  file = "analysis/psdan02_dm concordance pooled_or.csv",
  row.names = FALSE
)

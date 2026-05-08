# Purpose: Build spouse dyads per imputation and compute time-to-event variables.
# Steps: Split by sex, join into dyads, flag at-risk cohorts, define helper summaries.
# Output: Saves the spouse dyad dataset list for analysis scripts.
rm(list=ls());gc();source(".Rprofile")

library(tibble)

carrs_recode_list <- readRDS(paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre02_recoded dataset_list.RDS"))

dyads_list <- vector("list", length(carrs_recode_list))

for (imp_num in seq_along(carrs_recode_list)) {
  carrs_recode <- carrs_recode_list[[imp_num]]

  individual <- carrs_recode %>%
    dplyr::filter(biocohort == 1) %>%
    mutate(
      DateDMbiomarker = as.Date(DateDMbiomarker),
      visit_date0 = as.Date(visit_date0),
      time_to_dm = as.numeric(DateDMbiomarker - visit_date0) / 365.2
    )

  hh_biocohort <- individual %>%
    group_by(hhid) %>%
    summarise(
      sum_biocohort = sum(biocohort, na.rm = TRUE),
      count_hh = n(),
      .groups = "drop"
    )

  husb <- individual %>%
    dplyr::filter(sex == "male") %>%
    dplyr::distinct(hhid, .keep_all = TRUE) %>%
    dplyr::select(hhid, age, bmi, waist_cm, educcat, educcat2, employocccat, bmicat, bmicat2,
                  alc_overall, smk_overall, famhx_htn, famhx_dm, famhx_cvd, htn, dm, chd, cva, ckd,
                  multimorbiditycat, dm_biomarker0, event_DMbiomarker, DateDMbiomarker, time_to_dm, biocohort) %>%
    dplyr::rename_with(~ paste0(.x, "_husb"), -hhid)

  wife <- individual %>%
    dplyr::filter(sex == "female") %>%
    dplyr::distinct(hhid, .keep_all = TRUE) %>%
    dplyr::select(hhid, age, bmi, waist_cm, educcat, educcat2, employocccat, bmicat, bmicat2,
                  alc_overall, smk_overall, famhx_htn, famhx_dm, famhx_cvd, htn, dm, chd, cva, ckd,
                  multimorbiditycat, dm_biomarker0, event_DMbiomarker, DateDMbiomarker, time_to_dm, biocohort) %>%
    dplyr::rename_with(~ paste0(.x, "_wife"), -hhid)

  dyads0 <- husb %>% dplyr::inner_join(wife, by = "hhid")

  dyads_list[[imp_num]] <- dyads0 %>%
    mutate(
      dyad_at_risk = case_when(
        dm_biomarker0_husb == 1 & dm_biomarker0_wife == 1 ~ NA_real_,
        TRUE ~ 1
      ),
      female_at_risk = case_when(
        dm_biomarker0_wife == 1 ~ NA_real_,
        rowSums(is.na(cbind(age_wife, bmi_wife, educcat2_wife, famhx_dm_wife, multimorbiditycat_wife))) > 0 ~ NA_real_,
        TRUE ~ 1
      ),
      male_at_risk = case_when(
        dm_biomarker0_husb == 1 ~ NA_real_,
        rowSums(is.na(cbind(age_husb, bmi_husb, educcat2_husb, famhx_dm_husb, multimorbiditycat_husb,
                            smk_overall_husb, alc_overall_husb))) > 0 ~ NA_real_,
        TRUE ~ 1
      )
    )
}

saveRDS(
  dyads_list,
  paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS")
)

# Pooled spouse descriptives using Rubin's rules
mean_table <- function(df, var_name, group_var) {
  df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::group_modify(~ {
      x <- .x[[var_name]]
      nobs <- nrow(.x)
      nmiss <- sum(is.na(x))
      n_nonmiss <- nobs - nmiss
      mean_x <- if (n_nonmiss > 0) mean(x, na.rm = TRUE) else NA_real_
      sum_x <- if (n_nonmiss > 0) sum(x, na.rm = TRUE) else NA_real_

      if (n_nonmiss > 1) {
        sd_x <- stats::sd(x, na.rm = TRUE)
        se_mean <- sd_x / sqrt(n_nonmiss)
      } else {
        se_mean <- NA_real_
      }

      tibble::tibble(
        `Variable Name` = var_name,
        N = nobs,
        Mean = mean_x,
        Sum = sum_x,
        `Std Error of Mean` = se_mean,
        NMiss = nmiss
      )
    }) %>%
    dplyr::ungroup()
}

ratio_table <- function(df, num_var, denom_var, group_var) {
  df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::group_modify(~ {
      x <- .x[[num_var]]
      y <- .x[[denom_var]]
      nobs <- nrow(.x)
      complete <- stats::complete.cases(x, y)
      n_complete <- sum(complete)

      if (n_complete > 1) {
        x_c <- x[complete]
        y_c <- y[complete]
        mean_x <- mean(x_c)
        mean_y <- mean(y_c)
        ratio <- mean_x / mean_y
        z <- x_c - ratio * y_c
        var_z <- stats::var(z)
        se_ratio <- sqrt(var_z / (n_complete * mean_y^2))
      } else {
        ratio <- NA_real_
        se_ratio <- NA_real_
      }

      tibble::tibble(
        `Numerator Variable` = num_var,
        `Denominator Variable` = denom_var,
        N = nobs,
        Ratio = ratio,
        StdErr = se_ratio
      )
    }) %>%
    dplyr::ungroup()
}

spouse_descriptives_one <- function(dyads) {
  husb_events <- dyads %>%
    dplyr::filter(dm_biomarker0_husb == 0) %>%
    mean_table("event_DMbiomarker_husb", "dm_biomarker0_wife")

  husb_rates <- dyads %>%
    dplyr::filter(dm_biomarker0_husb == 0) %>%
    ratio_table("event_DMbiomarker_husb", "time_to_dm_husb", "dm_biomarker0_wife")

  wife_events <- dyads %>%
    dplyr::filter(dm_biomarker0_wife == 0) %>%
    mean_table("event_DMbiomarker_wife", "dm_biomarker0_husb")

  wife_rates <- dyads %>%
    dplyr::filter(dm_biomarker0_wife == 0) %>%
    ratio_table("event_DMbiomarker_wife", "time_to_dm_wife", "dm_biomarker0_husb")

  rates <- dplyr::bind_rows(wife_rates, husb_rates) %>%
    dplyr::mutate(spouse_diabetes = dplyr::coalesce(dm_biomarker0_husb, dm_biomarker0_wife))

  events <- dplyr::bind_rows(wife_events, husb_events) %>%
    dplyr::mutate(spouse_diabetes = dplyr::coalesce(dm_biomarker0_husb, dm_biomarker0_wife))

  events %>%
    dplyr::inner_join(
      rates,
      by = c("spouse_diabetes" = "spouse_diabetes", "Variable Name" = "Numerator Variable", "N")
    ) %>%
    mutate(
      population = case_when(
        `Variable Name` == "event_DMbiomarker_wife" ~ "Women",
        `Variable Name` == "event_DMbiomarker_husb" ~ "Men",
        TRUE ~ NA_character_
      ),
      dm_num = N,
      dm_events = Sum,
      dm_risk = Mean,
      dm_risk_se = `Std Error of Mean`,
      dm_incidence = Ratio,
      dm_incidence_se = StdErr
    ) %>%
    select(population, spouse_diabetes, dm_num, dm_events, dm_risk, dm_risk_se, dm_incidence, dm_incidence_se)
}

pool_rubin <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]
  se <- se[ok]
  m <- length(est)
  if (m == 0) {
    return(list(qbar = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
  }
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

valid_dyads <- dyads_list[!vapply(dyads_list, is.null, logical(1))]
spouse_all <- dplyr::bind_rows(lapply(valid_dyads, spouse_descriptives_one), .id = "imp")

get_pool_val <- function(x, name) {
  if (is.null(x) || is.null(x[[name]])) NA_real_ else x[[name]]
}

spouse_pooled <- spouse_all %>%
  dplyr::group_by(population, spouse_diabetes) %>%
  dplyr::reframe(
    dm_num = round(mean(dm_num, na.rm = TRUE)),
    dm_events = mean(dm_events, na.rm = TRUE),
    risk_pool = list(pool_rubin(dm_risk, dm_risk_se)),
    inc_pool = list(pool_rubin(dm_incidence, dm_incidence_se))
  ) %>%
  dplyr::mutate(
    dm_risk = vapply(risk_pool, get_pool_val, numeric(1), "qbar"),
    dm_incidence_100 = vapply(inc_pool, get_pool_val, numeric(1), "qbar") * 100,
    dm_incidence_100_lcl = vapply(inc_pool, get_pool_val, numeric(1), "lcl") * 100,
    dm_incidence_100_ucl = vapply(inc_pool, get_pool_val, numeric(1), "ucl") * 100
  ) %>%
  dplyr::select(
    population, spouse_diabetes, dm_num, dm_incidence_100, dm_incidence_100_ucl,
    dm_incidence_100_lcl, dm_events, dm_risk
  )

write.csv(
  spouse_pooled,
  paste0(path_spouses_diabetes_folder, "/Diabetes risk and incidence in spouse dyads.csv"),
  row.names = FALSE
)

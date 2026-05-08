# Purpose: Construct complete-case spouse dyads with time-to-event variables.
# Steps: Split by sex, join into dyads, flag at-risk cohorts.
# Output: Saves psdcpre03_spouse dyad dataset.RDS.
rm(list=ls());gc();source(".Rprofile")

carrs_recode <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre02_recoded dataset.RDS"))

individual <- carrs_recode %>% 
  dplyr::filter(biocohort == 1) %>% 
  mutate(time_to_dm = (as.numeric(DateDMbiomarker-visit_date0)/365.2))

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


dyads <- dyads0 %>%
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

saveRDS(dyads, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_spouse_dyad_dataset.RDS"))


# SAS-like mean and ratio outputs (SRS assumptions for SE/CI)
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
        tcrit <- stats::qt(0.975, df = n_nonmiss - 1)
        lower <- mean_x - tcrit * se_mean
        upper <- mean_x + tcrit * se_mean
      } else {
        se_mean <- NA_real_
        lower <- NA_real_
        upper <- NA_real_
      }

      se_sum <- if (!is.na(se_mean)) se_mean * n_nonmiss else NA_real_

      tibble::tibble(
        `Variable Name` = var_name,
        N = nobs,
        Mean = mean_x,
        Sum = sum_x,
        LowerCLMean = lower,
        `Std Error of Mean` = se_mean,
        `Std Error of Sum` = se_sum,
        UpperCLMean = upper,
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
        tcrit <- stats::qt(0.975, df = n_complete - 1)
        lower <- ratio - tcrit * se_ratio
        upper <- ratio + tcrit * se_ratio
      } else {
        ratio <- NA_real_
        se_ratio <- NA_real_
        lower <- NA_real_
        upper <- NA_real_
      }

      tibble::tibble(
        `Numerator Variable` = num_var,
        `Denominator Variable` = denom_var,
        N = nobs,
        Ratio = ratio,
        LowerCL = lower,
        StdErr = se_ratio,
        UpperCL = upper
      )
    }) %>%
    dplyr::ungroup()
}

# Event rates for husbands (overall + by wife's baseline status)
husb_events <- dyads %>%
  dplyr::filter(dm_biomarker0_husb == 0) %>%
  mean_table("event_DMbiomarker_husb", "dm_biomarker0_wife")

husb_rates <- dyads %>%
  dplyr::filter(dm_biomarker0_husb == 0) %>%
  ratio_table("event_DMbiomarker_husb", "time_to_dm_husb", "dm_biomarker0_wife")

# Event rates for wives (overall + by husband's baseline status)
wife_events <- dyads %>%
  dplyr::filter(dm_biomarker0_wife == 0) %>%
  mean_table("event_DMbiomarker_wife", "dm_biomarker0_husb")

wife_rates <- dyads %>%
  dplyr::filter(dm_biomarker0_wife == 0) %>%
  ratio_table("event_DMbiomarker_wife", "time_to_dm_wife", "dm_biomarker0_husb")

# Creating a single rates dataset
rates <- dplyr::bind_rows(wife_rates, husb_rates) %>%
  dplyr::mutate(spouse_diabetes = dplyr::coalesce(dm_biomarker0_husb, dm_biomarker0_wife))

# Creating a single events dataset
events <- dplyr::bind_rows(wife_events, husb_events) %>%
  dplyr::mutate(spouse_diabetes = dplyr::coalesce(dm_biomarker0_husb, dm_biomarker0_wife))

# Merge and derive descriptive fields
spouse_descriptives <- events %>%
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
    dm_incidence_100 = Ratio * 100,
    dm_incidence_100_ucl = UpperCL * 100,
    dm_incidence_100_lcl = LowerCL * 100,
    dm_num = N,
    dm_events = round(Sum, 5),
    dm_risk = Mean
  ) %>%
  select(
    population, spouse_diabetes, dm_num, dm_incidence_100, dm_incidence_100_ucl,
    dm_incidence_100_lcl, dm_events, dm_risk
  )

write.csv(
  spouse_descriptives,
  paste0(path_spouses_diabetes_folder, "/cca results/Diabetes risk and incidence in spouse dyads cca.csv"),
  row.names = FALSE
)

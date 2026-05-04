rm(list=ls());gc();source(".Rprofile")

library(dplyr)

# Table 3 pooled across imputed dyads

dyads_list <- readRDS(paste0(path_spouses_diabetes_folder,"/working/preprocessing/psdpre03_spouse dyad dataset_list.RDS"))
dyads_list <- dyads_list[!vapply(dyads_list, is.null, logical(1))]

pool_rubin <- function(est, se) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]
  se <- se[ok]
  m <- length(est)
  if (m == 0) return(list(qbar = NA_real_, se = NA_real_))
  qbar <- mean(est)
  ubar <- mean(se^2)
  b <- stats::var(est)
  if (is.na(b)) b <- 0
  tvar <- ubar + (1 + 1 / m) * b
  list(qbar = qbar, se = sqrt(tvar))
}

one_imp <- function(dyads) {
  women_df <- dyads %>%
    dplyr::filter(dm_biomarker0_wife == 0) %>%
    dplyr::mutate(
      partner_status = dplyr::case_when(
        dm_biomarker0_husb == 0 ~ "No diabetes",
        dm_biomarker0_husb == 1 ~ "Diabetes",
        TRUE ~ NA_character_
      ),
      outcome_status = dplyr::case_when(
        event_DMbiomarker_wife == 0 ~ "No diabetes",
        event_DMbiomarker_wife == 1 ~ "Diabetes",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

  men_df <- dyads %>%
    dplyr::filter(dm_biomarker0_husb == 0) %>%
    dplyr::mutate(
      partner_status = dplyr::case_when(
        dm_biomarker0_wife == 0 ~ "No diabetes",
        dm_biomarker0_wife == 1 ~ "Diabetes",
        TRUE ~ NA_character_
      ),
      outcome_status = dplyr::case_when(
        event_DMbiomarker_husb == 0 ~ "No diabetes",
        event_DMbiomarker_husb == 1 ~ "Diabetes",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(partner_status), !is.na(outcome_status))

  build <- function(df, sex_label) {
    df %>%
      dplyr::count(partner_status, outcome_status, name = "n") %>%
      dplyr::group_by(partner_status) %>%
      dplyr::mutate(
        partner_n = sum(n),
        p = n / partner_n,
        se = sqrt(p * (1 - p) / partner_n)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sex = sex_label)
  }

  dplyr::bind_rows(build(women_df, "Women"), build(men_df, "Men"))
}

all_imp <- dplyr::bind_rows(lapply(dyads_list, one_imp), .id = "imp")

pooled <- all_imp %>%
  dplyr::group_by(sex, partner_status, outcome_status) %>%
  dplyr::reframe(
    n_mean = mean(n, na.rm = TRUE),
    partner_n = mean(partner_n, na.rm = TRUE),
    pool = list(pool_rubin(p, se))
  ) %>%
  dplyr::mutate(
    p = vapply(pool, `[[`, numeric(1), "qbar")
  )

women_tbl <- pooled %>%
  dplyr::filter(sex == "Women") %>%
  dplyr::mutate(
    baseline_group = paste0("Women without diabetes at baseline (n=", round(mean(partner_n, na.rm = TRUE)), ")"),
    partner_group = paste0(partner_status, " (n=", round(partner_n), ")")
  ) %>%
  dplyr::select(baseline_group, partner_group, outcome_status, p) %>%
  tidyr::pivot_wider(names_from = outcome_status, values_from = p)

men_tbl <- pooled %>%
  dplyr::filter(sex == "Men") %>%
  dplyr::mutate(
    baseline_group = paste0("Men without diabetes at baseline (n=", round(mean(partner_n, na.rm = TRUE)), ")"),
    partner_group = paste0(partner_status, " (n=", round(partner_n), ")")
  ) %>%
  dplyr::select(baseline_group, partner_group, outcome_status, p) %>%
  tidyr::pivot_wider(names_from = outcome_status, values_from = p)

final_tbl <- dplyr::bind_rows(women_tbl, men_tbl) %>%
  dplyr::rename(
    `Individual's status at baseline` = baseline_group,
    `Individuals with following partner's status at baseline` = partner_group,
    `No diabetes, %` = `No diabetes`,
    `Diabetes, %` = Diabetes
  ) %>%
  dplyr::mutate(
    `No diabetes, %` = sprintf("%.1f", 100 * `No diabetes, %`),
    `Diabetes, %` = sprintf("%.1f", 100 * `Diabetes, %`)
  )

utils::write.csv(
  final_tbl,
  file = "analysis/psdan03_followup_risk_by_spouse_baseline.csv",
  row.names = FALSE
)

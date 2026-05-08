# Purpose: Summarize person-time by sex, event status, and spouse baseline DM exposure.
# Steps: Filter at-risk dyads, group by outcome/exposure, summarize time_to_dm.
# Output: Writes a person-time summary table to cca/analysis/.
rm(list=ls());gc();source(".Rprofile")

# Load complete-case dyad dataset

dyads <- readRDS(
  paste0(path_spouses_diabetes_folder, "/working/cca/preprocessing/psdcpre03_spouse_dyad_dataset.RDS")
)

summarize_person_time <- function(df, sex_label, time_var, event_var, spouse_var, baseline_var) {
  df %>%
    dplyr::filter(
      dyad_at_risk == 1,
      .data[[baseline_var]] == 0,
      !is.na(.data[[time_var]]),
      !is.na(.data[[event_var]]),
      !is.na(.data[[spouse_var]])
    ) %>%
    dplyr::mutate(
      event_status = .data[[event_var]],
      exposure_status = .data[[spouse_var]]
    ) %>%
    dplyr::group_by(event_status, exposure_status) %>%
    dplyr::summarise(
      sum_person_time = sum(.data[[time_var]], na.rm = TRUE),
      mean_person_time = mean(.data[[time_var]], na.rm = TRUE),
      min_person_time = min(.data[[time_var]], na.rm = TRUE),
      max_person_time = max(.data[[time_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(`Target population` = sex_label) %>%
    dplyr::select(
      `Target population`,
      `Event status at endline` = event_status,
      `Exposure status (sp has dm or not)` = exposure_status,
      `Sum of all person time` = sum_person_time,
      `Mean person time` = mean_person_time,
      `Min person time` = min_person_time,
      `Max person time` = max_person_time
    )
}

women_tbl <- summarize_person_time(
  dyads,
  sex_label = "Women",
  time_var = "time_to_dm_wife",
  event_var = "event_DMbiomarker_wife",
  spouse_var = "dm_biomarker0_husb",
  baseline_var = "dm_biomarker0_wife"
)

men_tbl <- summarize_person_time(
  dyads,
  sex_label = "Men",
  time_var = "time_to_dm_husb",
  event_var = "event_DMbiomarker_husb",
  spouse_var = "dm_biomarker0_wife",
  baseline_var = "dm_biomarker0_husb"
)

person_time_tbl <- dplyr::bind_rows(women_tbl, men_tbl)

utils::write.csv(
  person_time_tbl,
  file = "cca/analysis/psdcan07_person_time_table.csv",
  row.names = FALSE
)

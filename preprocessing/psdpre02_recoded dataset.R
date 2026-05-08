# Purpose: Create analysis flags and recoded covariates on imputed wide datasets.
# Steps: Identify valid dyads/follow-up, derive incident DM timing, recode covariates.
# Output: Saves a recoded dataset list for dyad construction.
rm(list=ls());gc();source(".Rprofile")

carrs_wide_list <- readRDS(paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre01_wide analytic dataset_list.RDS"))

recoded_list <- vector("list", length(carrs_wide_list))

for (imp_num in seq_along(carrs_wide_list)) {
  df <- carrs_wide_list[[imp_num]]

  # Household-level summary
  hh_summary <- df %>%
    dplyr::filter(!is.na(hhid)) %>%
    group_by(hhid) %>%
    summarise(
      sum_spouse = sum(spousedyad_new, na.rm = TRUE),
      count_hh = n(),
      .groups = "drop"
    )

  # Individual-level: create inclusion/exclusion indicators
  carrs_wide1 <- df %>%
    left_join(hh_summary, by = "hhid") %>%
    mutate(
      valid_dyad = case_when(
        sum_spouse == 2 & count_hh == 2 ~ 1,
        TRUE ~ 0
      ),
      baseline_bio_all = case_when(
        dm_biomarker_nomiss0 %in% c(0, 1) ~ 1,
        TRUE ~ 0
      ),
      baseline_bio_one = case_when(
        !is.na(fpg0) | !is.na(hba1c0) ~ 1,
        TRUE ~ 0
      ),
      baseline_self = case_when(
        !is.na(dm) ~ 1,
        TRUE ~ 0
      ),
      nfollowup = rowSums(cbind(visitcomplete1, visitcomplete2, visitcomplete3), na.rm = TRUE),
      onefollowup = case_when(
        nfollowup >= 1 ~ 1,
        TRUE ~ 0
      ),
      completefollowup = case_when(
        carrs == 1 & nfollowup == 3 ~ 1,
        carrs == 2 & visitcomplete1 == 1 ~ 1,
        TRUE ~ 0
      ),
      onefollowupbio = case_when(
        carrs == 1 & rowSums(is.na(cbind(dm_biomarker1, dm_biomarker2, dm_biomarker3))) < 3 ~ 1,
        carrs == 2 & !is.na(dm_biomarker1) ~ 1,
        TRUE ~ 0
      ),
      completefollowupbio = case_when(
        carrs == 1 & rowSums(is.na(cbind(dm_biomarker1, dm_biomarker2, dm_biomarker3))) == 0 ~ 1,
        carrs == 2 & !is.na(dm_biomarker1) ~ 1,
        TRUE ~ 0
      ),
      biocohort = case_when(
        valid_dyad == 1 & baseline_bio_one == 1 & onefollowupbio == 1 ~ 1,
        TRUE ~ 0
      ),
      valid_at_risk = case_when(
        valid_dyad == 1 & dm_biomarker0 == 0 & onefollowupbio == 1 ~ 1,
        TRUE ~ 0
      )
    )

  # Creating event dates in the wide dataset (vectorized)
  follow_mat <- as.matrix(carrs_wide1[, c("dm_biomarker1", "dm_biomarker2", "dm_biomarker3")])
  visit_mat <- as.matrix(carrs_wide1[, c("visit_date1", "visit_date2", "visit_date3")])

  any_followup <- rowSums(!is.na(follow_mat)) > 0
  any_event <- rowSums(follow_mat == 1, na.rm = TRUE) > 0

  first_idx <- ifelse(any_event, max.col(follow_mat == 1, ties.method = "first"), NA_integer_)
  last_idx <- ifelse(any_followup, max.col(!is.na(visit_mat), ties.method = "last"), NA_integer_)

  event_DMbiomarker <- ifelse(
    carrs_wide1$dm_biomarker0 == 0,
    ifelse(!any_followup, NA_real_, ifelse(any_event, 1, 0)),
    NA_real_
  )

  first_dates <- as.Date(visit_mat[cbind(seq_len(nrow(carrs_wide1)), first_idx)], origin = "1970-01-01")
  last_dates <- as.Date(visit_mat[cbind(seq_len(nrow(carrs_wide1)), last_idx)], origin = "1970-01-01")

  DateDMbiomarker <- ifelse(
    carrs_wide1$dm_biomarker0 == 0,
    ifelse(event_DMbiomarker == 1, first_dates,
           ifelse(event_DMbiomarker == 0, last_dates, as.Date(NA))),
    as.Date(NA)
  )

  incident_outcomes <- carrs_wide1 %>%
    mutate(
      firstvisit_DMbiomarker = first_idx,
      event_DMbiomarker = event_DMbiomarker,
      DateDMbiomarker = DateDMbiomarker
    )

  recoded_list[[imp_num]] <- incident_outcomes %>%
    mutate(
      educcat = case_when(
        educstat %in% c(7, 6, 5) ~ 1,
        educstat %in% c(3, 4) ~ 2,
        educstat %in% c(1, 2) ~ 3,
        TRUE ~ NA_real_
      ),
      # education category (binary)
      educcat2 = case_when(
        educstat %in% c(3, 4, 5, 6, 7) ~ 1,
        educstat %in% c(1, 2) ~ 2,
        TRUE ~ NA_real_
      ),
      employocccat = case_when(
        employ %in% c(2, 3)                        ~ "Not in the labor force, student/housewives",
        employ == 4                                ~ "Not in the labor force, retired",
        employ == 5                                ~ "Unemployed",
        employ == 1                                ~ "Employed",
        TRUE                                       ~ "Others"
      ),
      bmicat = case_when(
        is.na(bmi) ~ NA_real_,
        bmi < 25 ~ 0,
        bmi >= 25 & bmi < 30 ~ 1,
        bmi >= 30 ~ 2,
        TRUE ~ NA_real_
      ),
      bmicat2 = case_when(
        is.na(bmi) ~ NA_real_,
        bmi < 30 ~ 1,
        bmi >= 30 ~ 2,
        TRUE ~ NA_real_
      ),
      htn = replace_na(htn, 0),
      dm = replace_na(dm, 0),
      chd = replace_na(chd, 0),
      cva = replace_na(cva, 0),
      ckd = replace_na(ckd, 0),
      multimorbidity = htn + dm + chd + cva + ckd,
      multimorbiditycat = case_when(
        multimorbidity == 0 ~ 0,
        multimorbidity == 1 ~ 1,
        multimorbidity >= 2 ~ 2,
        TRUE ~ NA_real_
      ),
      multimorbiditycat2 = case_when(
        multimorbiditycat == 0 ~ 0,
        multimorbiditycat > 0 ~ 1,
        TRUE ~ NA_real_
      )
    )
  }

saveRDS(
  recoded_list,
  paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre02_recoded dataset_list.RDS")
)


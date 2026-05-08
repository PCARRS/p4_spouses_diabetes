# Purpose: Apply inclusion flags and derive incident DM timing for complete-case data.
# Steps: Identify valid dyads/follow-up, compute event dates, recode covariates.
# Output: Saves the complete-case recoded dataset for dyad construction.
rm(list=ls());gc();source(".Rprofile")

carrs_wide <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre01_wide analytic dataset.RDS"))

# Household-level summary
hh_summary <- carrs_wide %>%
	dplyr::filter(!is.na(hhid)) %>%
	group_by(hhid) %>%
	summarise(
		sum_spouse = sum(spousedyad_new, na.rm = TRUE),
		count_hh = n(),
		.groups = "drop"
	)

# Individual-level: create inclusion/exclusion indicators
carrs_wide1 <- carrs_wide %>%
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

# Creating event dates in the wide dataset
incident_outcomes <- carrs_wide1 %>%
	rowwise() %>%
	mutate(
		# First follow-up visit (1-3) where DM biomarker indicates incident DM
		firstvisit_DMbiomarker = case_when(
			dm_biomarker0 == 0 ~ {
				dm_followup <- c_across(dm_biomarker1:dm_biomarker3)
				idx <- which(dm_followup == 1)
				if (length(idx) == 0) NA_real_ else idx[1]
			},
			TRUE ~ NA_real_
		),
		# Event status: 1 if incident DM occurs, 0 if not, NA if no follow-up data
		event_DMbiomarker = case_when(
			dm_biomarker0 == 0 ~ {
				dm_followup <- c_across(dm_biomarker1:dm_biomarker3)
				if (all(is.na(dm_followup))) {
					NA_real_
				} else if (!is.na(firstvisit_DMbiomarker)) {
					1
				} else {
					0
				}
			},
			TRUE ~ NA_real_
		),
		# Event date: visit date at first incident DM; otherwise last non-missing visit date
		DateDMbiomarker = case_when(
			dm_biomarker0 == 0 ~ {
				visit_dates <- c_across(visit_date1:visit_date3)
				if (is.na(event_DMbiomarker)) {
					as.Date(NA)
				} else if (event_DMbiomarker == 1) {
					visit_dates[firstvisit_DMbiomarker]
				} else {
					last_idx <- which(!is.na(visit_dates))
					if (length(last_idx) == 0) as.Date(NA) else visit_dates[max(last_idx)]
				}
			},
			TRUE ~ as.Date(NA)
		)
	) %>%
	ungroup()


carrs_recode <- incident_outcomes %>%
  mutate(
    # education category
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
    # employment & occupation category
    employocccat = case_when(
      employ == 1 & occ %in% c(3, 4, 5) ~ 4,
      employ == 1 & occ %in% c(1, 2) ~ 5,
      employ %in% c(2, 3) ~ 1,
      employ == 4 ~ 2,
      employ == 5 ~ 3,
      TRUE ~ NA_real_
    ),
    # BMI categories
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
    # multimorbidity (treat missing indicators as 0 for counting)
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

saveRDS(carrs_recode, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre02_recoded_dataset.RDS"))


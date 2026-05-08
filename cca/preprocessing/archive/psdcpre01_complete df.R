# Archive: Build complete-case spouse dyads from long data with imputed dates.
# Steps: Trim non-response visits, impute dates, enforce one male/one female per household.
# Output: Saves the long spouse complete-case dataset.
rm(list=ls());gc();source(".Rprofile")

# Input harmonized CARRS dataset (Exclude - Karachi site), N = 21,862
carrs_df <- readRDS(paste0(path_spouses_bmi_change_folder,"/working/preprocessing/psbpre02_carrs harmonized data.RDS")) %>% 
  # reason: indication of why participant did not complete interview.
  # if reason for not participating is missing, we assume the visit was completed;
  mutate(
    visitcomplete = case_when(
      is.na(reason) ~ 1,
      TRUE ~ 0
    )) %>% 
  mutate(
    visit_date = case_when(
      visitcomplete == 1 ~ doi,
      TRUE ~ as.Date(NA)
    )) %>% 
  # define diabetes event
  mutate(
    dm_biomarker = case_when(
      visitcomplete == 1 & is.na(fpg) & is.na(hba1c) ~ NA_real_,
      visitcomplete == 1 & (fpg >= 126 | hba1c >= 6.5 | dm == 1 | dm_allo == 1) ~ 1,
      visitcomplete == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    dm_biomarker_nomiss = case_when(
      visitcomplete == 1 & (is.na(fpg) | is.na(hba1c)) ~ NA_real_,
      visitcomplete == 1 & (fpg >= 126 | hba1c >= 6.5 | dm == 1 | dm_allo == 1) ~ 1,
      visitcomplete == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    dm_self = case_when(
      visitcomplete == 1 & (dm == 1 | dm_allo == 1) ~ 1,
      visitcomplete == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    # count missing DM biomarkers
    missingbiomarkern = case_when(
      visitcomplete == 1 ~ rowSums(is.na(cbind(fpg, hba1c))),
      TRUE ~ NA_real_
    )
  ) %>% 
  # valid dyad: both members in the same household are flagged spouseyad_new == 1
  group_by(carrs, hhid) %>%
  mutate(valid_dyad = case_when(
    n_distinct(pid[spousedyad_new == 1], na.rm = TRUE) == 2 &
      n_distinct(pid, na.rm = TRUE) == 2 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup()


# Step 0: keep main visits only

carrs_main <- carrs_df %>% 
  dplyr::filter(!is.na(visit)) %>%
  group_by(carrs, pid) %>%
  mutate(nfollowup = n_distinct(visit, na.rm = TRUE)) %>%
  mutate(
    complete_followup = case_when(
      carrs == 1 & nfollowup == 4 ~ 1,
      carrs == 2 & nfollowup == 2 ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  ungroup()


# Step 1: Keep valid spouses only

spouse_df <- carrs_main %>% 
  # N = 13,207, 6,712 spouses
  dplyr::filter(spousedyad_new == 1) %>% 
  # valid dyad, N = 12,990
  dplyr::filter(valid_dyad == 1)


# Step 2: Biomarker availability

# Has at least one biomarker for diabetes at baseline, N = 11,987
dmbio_either <- spouse_df %>% 
  dplyr::filter(fup == 0, missing_biomarker %in% c(0,1))

# Diabetes biomarker available at least one time point after baseline, N = 9,058
dmbio_fup_ids <- spouse_df %>% 
  dplyr::filter(pid %in% dmbio_either$pid, fup > 0) %>%
  group_by(carrs, pid) %>%
  summarise(has_fup_dm_biomarker = any(!is.na(dm_biomarker)), .groups = "drop") %>%
  dplyr::filter(has_fup_dm_biomarker) %>%
  select(carrs, pid)

dmbio_fup <- spouse_df %>% 
  semi_join(dmbio_fup_ids, by = c("carrs", "pid"))

# Check number of pid per hhid
dmbio_fup %>%
  distinct(hhid, pid) %>%
  count(hhid, name = "n_pid") %>%
  count(n_pid, name = "Freq")



# Recode key covariates
recoded_df <- dmbio_fup %>%   
  arrange(hhid, pid, carrs, fup) %>%
  # Nielsen et al. 2023 (based on the mean of the second two of three blood pressure measurements)
  mutate(
    sbp = rowMeans(select(., sbp2, sbp3), na.rm = TRUE),
    dbp = rowMeans(select(., dbp2, dbp3), na.rm = TRUE)
  ) %>%
  # define disease indicators
  mutate(
    overweight = case_when(
      is.na(bmi) ~ NA_real_,
      bmi >= 25 ~ 1,
      TRUE ~ 0
    ),
    hypertension = case_when(
      is.na(sbp) & is.na(dbp) & is.na(htn) & is.na(htn_med) & is.na(htn_allo) ~ NA_real_,
      sbp > 140 | dbp > 90 | htn == 1 | htn_med == 1 | htn_allo == 1 ~ 1,
      TRUE ~ 0
    ),
    high_tg = case_when(
      is.na(tg) ~ NA_real_,
      tg > 150 ~ 1,
      TRUE ~ 0
    ),
    chd = case_when(
      is.na(chd_allo) ~ NA_real_,
      chd_allo == 1 ~ 1,
      TRUE ~ chd
    )
  ) %>% 
  mutate(morbidity_number = rowSums(across(c(hypertension, dm_biomarker, chd, cva, ckd), ~ .x == 1), na.rm = TRUE)) %>% 
  mutate(
    morbidity_category = case_when(
      morbidity_number == 0 ~ "None",
      morbidity_number == 1 ~ "Single morbidity",
      TRUE                  ~ "Multimorbidity"
    )) %>% 
  # Unknown/NA to 0 
  mutate(famhx_cvd = case_when(is.na(famhx_cvd) ~ 0, 
                               TRUE ~ famhx_cvd),
         famhx_htn = case_when(is.na(famhx_htn) ~ 0, 
                               TRUE ~ famhx_htn),
         famhx_dm = case_when(is.na(famhx_dm) ~ 0, 
                              TRUE ~ famhx_dm)) %>% 
  mutate(
    alc_overall = case_when(
      alc_curr == 1 ~ 1,  # current (follow-up)
      alc_curr == 0 ~ 3,  # assume never
      alc_often %in% c(1, 2) | alc_overall == 1 ~ 1,  # current (baseline)
      alc_often %in% c(3, 4) | alc_overall == 2 ~ 2,  # former
      is.na(alc_often) | alc_overall == 0 ~ 0,        # never
      TRUE ~ alc_overall
    ),
    smk_overall = case_when(
      smk_curr == 1 | smk_smoke_freq %in% c(1, 2) | smk_chew_freq %in% c(1, 2) ~ 1, # current
      smk_curr == 0 | smk_ever == 0 | smk_overall == 0 ~ 0,                         # never
      smk_ever == 1 | smk_overall == 2 ~ 2,                                        # former
      TRUE ~ smk_overall
    )
  ) %>% 
  mutate(
    bmi_category = case_when(
      bmi >= 15 & bmi < 25 ~ "Underweight or normal weight",
      bmi >= 25 & bmi < 30 ~ "Overweight",
      bmi >= 30 ~ "Obese",
      TRUE ~ NA_character_
    )
  ) %>%
  # Define DM biomarker event for survival analysis
  group_by(carrs, pid) %>%
  mutate(
    event_dm_biomarker = case_when(
      all(is.na(dm_biomarker[fup > 0])) ~ NA_real_,
      any(dm_biomarker == 1 & fup > 0, na.rm = TRUE) ~ 1,
      TRUE ~ 0
    ),
    # Date of first DM biomarker event (or NA if no event)
    date_dm_biomarker = if_else(
      event_dm_biomarker == 1,
      min(visit_date[dm_biomarker == 1 & fup > 0], na.rm = TRUE),
      as.Date(NA)
    ),
    # Baseline DM biomarker status (fup == 0)
    dm_biomarker0 = case_when(
      fup == 0 ~ dm_biomarker,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()


# Step 4: 

# one row per person

# transpose to wide format, one row per spouse

value_cols <- setdiff(
  names(recoded_df),
  c("hhid", "sex", "carrs", "fup")    # drop your id‐cols here
)
# unique hhid: 1,035
dmbs_free_wide <- recoded_df %>%
  pivot_wider(
    id_cols    = c(hhid, carrs, fup),   
    names_from = sex,
    values_from = all_of(value_cols),
    names_glue = "{sex}_{.value}"
  ) %>% 
  # Keep only household-wave records where both spouses are present after widening.
  dplyr::filter(!is.na(female_pid), !is.na(male_pid))

baseline_wide <- dmbs_free_wide %>% 
  dplyr::filter(fup == 0) %>% 
  mutate(
    # Removing couples where both partners have diabetes at baseline
    dyad_at_risk = case_when(
      female_dm_biomarker == 1 & male_dm_biomarker == 1 ~ NA_real_,
      TRUE ~ 1
    ),
    # Female at risk: not diabetic at baseline AND complete covariates
    female_at_risk = case_when(
      female_dm_biomarker == 1 ~ NA_real_,
      if_any(
        all_of(paste0("female_", c("age", "bmi", "edu_category", "famhx_dm", "morbidity_category"))),
        is.na
      ) ~ NA_real_,
      TRUE ~ 1
    ),
    # Male at risk: not diabetic at baseline AND complete covariates
    male_at_risk = case_when(
      male_dm_biomarker == 1 ~ NA_real_,
      if_any(
        all_of(paste0("male_", c("age", "bmi", "edu_category", "famhx_dm", "morbidity_category", "smk_overall", "alc_overall"))),
        is.na
      ) ~ NA_real_,
      TRUE ~ 1
    )
  )



# Step 3: Free of diabetes at baseline (based on at least one biomarker, FPG or HbA1c), N = 6,673

female_dmfree <- dmbs_free_wide %>% 
  dplyr::filter(female_dm_biomarker0 == 0) %>% 
  distinct(carrs, hhid, female_pid, male_pid)

male_dmfree <- dmbs_free_wide %>% 
  dplyr::filter(male_dm_biomarker0 == 0) %>% 
  distinct(carrs, hhid, female_pid, male_pid)

ids <- bind_rows(female_dmfree,male_dmfree) %>% 
  distinct(carrs, hhid, female_pid, male_pid)

df <- baseline_wide %>% 
  semi_join(ids, by = c("carrs", "hhid", "female_pid", "male_pid"))

# Step 6: All other baseline covariates available

# N = 4,903
complete_pids <- dmbio_fup %>% 
  dplyr::filter(fup == 0) %>% 
  # delete rows with any of these missing
  dplyr::filter(
    if_all(
      all_of(c(
        "age", "sex", "edu_category", "employ_category", "hhincome",
        "bmi", "famhx_dm", "smk_overall"
      )),
      ~ !is.na(.)
    )
  ) %>% 
  distinct(pid)

complete_df <- dmbio_fup %>% 
  dplyr::filter(pid %in% complete_pids$pid)


saveRDS(dmbs_free, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre01_complete df.RDS"))


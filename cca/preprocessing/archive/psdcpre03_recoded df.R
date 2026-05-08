# Archive: Recode covariates for spouse dataset with imputed interview dates.
# Steps: Derive morbidity, smoking/alcohol, and BMI categories.
# Output: Saves psdcpre03_recoded_df.RDS.
rm(list=ls());gc();source(".Rprofile")

spouse_df <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre02_spouses with imputed date of interview.RDS"))

# Recode key covariates
recoded_df <- spouse_df %>%   
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
  )

saveRDS(recoded_df, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_recoded_df.RDS"))


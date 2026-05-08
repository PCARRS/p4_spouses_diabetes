# Purpose: Build complete-case analytic dataset (wide) from harmonized CARRS data.
# Steps: Define visit completion, derive diabetes indicators, merge baseline and spouse dyads.
# Output: Saves the complete-case wide analytic dataset for CCA.
rm(list=ls());gc();source(".Rprofile")

library(haven)

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
  ) 

# check percentage of missing in FPG and HbA1c
na_by_visit_carrs <- carrs_df %>%
  group_by(carrs, visit) %>%
  summarise(
    n_total = n(),
    nmiss_fpg = sum(is.na(fpg)),
    pctmiss_fpg = 100 * nmiss_fpg / n_total,
    nmiss_hba1c = sum(is.na(hba1c)),
    pctmiss_hba1c = 100 * nmiss_hba1c / n_total,
    nmiss_both = sum(is.na(fpg) & is.na(hba1c)),
    pctmiss_both = 100 * nmiss_both / n_total,
    nmiss_either = sum(is.na(fpg) | is.na(hba1c)),
    pctmiss_either = 100 * nmiss_either / n_total,
    .groups = "drop"
  ) %>%
  arrange(carrs, visit)

# A long dataset with only follow ups with biomarkers
carrs_df2 <- carrs_df %>% 
  dplyr::filter(!is.na(visit))

# Creating wide data sets for each of the key variables for diabetes incidence, N = 21,862
wide_outcomes <- carrs_df2 %>%
  select(
    pid, visit, fpg, hba1c, dm, dm_allo, dm_biomarker, dm_biomarker_nomiss,
    dm_self, visit_date, visitcomplete, missingbiomarkern
  ) %>%
  tidyr::pivot_wider(
    names_from = visit,
    values_from = c(
      fpg, hba1c, dm, dm_allo, dm_biomarker, dm_biomarker_nomiss,
      dm_self, visit_date, visitcomplete, missingbiomarkern
    ),
    names_glue = "{.value}{visit}"
  )

# N = 21,862
baseline <- read_sas(paste0(path_spouses_bmi_change_folder,"/working/raw/baseline_2025_0312.sas7bdat")) %>% 
  rename(carrs = CARRS) %>% 
  dplyr::select(
    # ID
    carrs,hhid,pid,
    # Demographic
    doi,dob,ceb,age,sex,site,educstat,employ,occ,hhincome,religion,caste,
    # Tobacco
    smk_ever,smk_curr,smk_overall,smk_exp,
    # Alcohol
    alc_often,alc_overall,
    # CVD
    htn, # HTN
    dm,dm_med,dm_allo,dm_rec, # DM
    hld, # heart disease
    chd,
    cva, # stroke
    ckd, # CKD
    cancer, # cancer
    # Family history
    famhx_htn,famhx_cvd,famhx_dm,
    # ANTHRO
    sbp1,sbp2,sbp3,dbp1,dbp2,dbp3,height_cm,weight_kg,bmi,waist_cm,hip_cm
    # Lab
    # fpg,tg,hba1c
  ) %>% 
  # keep Delhi, Chennai: 21,862
  dplyr::filter(site %in% c(1, 2)) %>% 
  mutate(
    fup = 0,
    site = case_when(site == 1 ~ "Chennai", TRUE ~ "Delhi"),
    sex = case_when(sex == 1 ~ "male", TRUE ~ "female"),
    bmi = case_when(is.na(bmi) ~ weight_kg / ((height_cm / 100) ^ 2), TRUE ~ bmi),
    year = as.integer(format(doi, "%Y"))
  )
# N = 21,862
spousedyads_clean <- readRDS(paste0(path_spouses_bmi_change_folder,"/working/preprocessing/spouseyads cleaned.RDS"))

carrs_wide <- baseline %>%
  left_join(wide_outcomes, by = "pid") %>%
  left_join(spousedyads_clean %>% 
            select(pid,hhid,spousedyad_new),
            by = c("pid","hhid"))

saveRDS(carrs_wide, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre01_wide analytic dataset.RDS"))


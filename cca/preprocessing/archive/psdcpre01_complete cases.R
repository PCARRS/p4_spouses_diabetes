# Archive: Define complete cases and identify eligible spouse dyads.
# Steps: Define responders, keep biomarker visits, impute dates, filter baseline covariates.
# Output: Saves complete-case spouse datasets for downstream steps.
rm(list=ls());gc();source(".Rprofile")

# ============================================================================
# Purpose: define complete cases and identify eligible spouse dyads for analysis.
# ============================================================================

carrs_df <- readRDS(paste0(path_spouses_bmi_change_folder,"/working/preprocessing/psbpre02_carrs harmonized data.RDS"))

# Step 1: Define responders - missing dropout reason

respond_df <- carrs_df %>%
  mutate(
    respond = case_when(
      (carrs == 1 & fup %in% c(1, 3, 5, 6)) | (carrs == 2 & fup == 1) ~ 0, # visits without biomarker measurements
      is.na(reason) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  arrange(carrs, pid, fup) %>%
  group_by(carrs, pid) %>%
  # Keep visits up to the last responded visit; trim only trailing non-response visits.
  mutate(
    visit_order = row_number(),
    last_respond_order = case_when(
      any(respond == 1) ~ max(visit_order[respond == 1]),
      TRUE ~ NA_integer_
    )
  ) %>%
  dplyr::filter(!is.na(last_respond_order), visit_order <= last_respond_order) %>%
  ungroup() %>%
  select(-visit_order, -last_respond_order)


# Step 2: Keep visits with biomarkers; define diabetes events

main_visits <- respond_df %>%
  dplyr::filter(
    (carrs == 1 & fup %in% c(0, 2, 4, 7)) |
      (carrs == 2 & fup %in% c(0, 2))
  ) %>% 
  mutate(
    dm_self_report = case_when(
      is.na(fpg) & is.na(hba1c) & is.na(dm) & is.na(dm_allo) ~ NA_real_,
      (dm == 1 | dm_allo == 1) & is.na(fpg) & is.na(hba1c) ~ 1,
      TRUE ~ 0
    ),
    dm_biomarker = case_when(
      respond == 0 ~ NA_real_,
      is.na(fpg) & is.na(hba1c) & is.na(dm) & is.na(dm_allo) ~ NA_real_,
      fpg >= 126 | hba1c >= 6.5 | dm == 1 | dm_allo == 1 ~ 1,
      TRUE ~ 0
    )
  )


# Step 7: Impute missing interview date with visit-specific midpoint
impute_date_df <- main_visits %>% 
  group_by(carrs, fup) %>%
  mutate(
    # Build visit-wave midpoint from observed DOI among responder visits only.
    doi_midpoint = {
      valid_doi <- doi[respond == 1 & !is.na(doi)]
      if (length(valid_doi) == 0) {
        doi[NA_integer_]
      } else {
        first_doi <- min(valid_doi, na.rm = TRUE)
        last_doi <- max(valid_doi, na.rm = TRUE)
        first_doi + as.numeric(last_doi - first_doi) / 2
      }
    },
    # Keep DOI missing for non-responders; impute only responder visits with missing DOI.
    doi = case_when(
      respond != 1 ~ as.Date(NA),
      respond == 1 & is.na(doi) ~ doi_midpoint,
      TRUE ~ doi
    )
  ) %>%
  ungroup() %>%
  select(-doi_midpoint) 

# Define inter-visit mid-point date - midpoint of the present visit and the last visit
inter_visit_df <- impute_date_df %>% 
  arrange(carrs, pid, fup) %>%
  group_by(carrs, pid) %>%
  mutate(
    baseline_date = case_when(
      any(!is.na(doi)) ~ first(doi[!is.na(doi)]),
      TRUE ~ as.Date(NA)
    ),
    doi_last_obs = doi
  ) %>%
  tidyr::fill(doi_last_obs, .direction = "down") %>%
  mutate(
    prior_visit_date = dplyr::lag(doi_last_obs),
    inter_visit_midpoint = case_when(
      !is.na(doi) & !is.na(prior_visit_date) ~
        prior_visit_date + as.numeric(doi - prior_visit_date) / 2,
      TRUE ~ as.Date(NA)
    )
  ) %>%
  select(-doi_last_obs) %>%
  ungroup()


# Step 1: Keep spouses only, N = 13,207, 6,712 spouses
spouse_df <- inter_visit_df %>% 
  dplyr::filter(spousedyad_new == 1) 

# Step 2: Exclude missing baseline covariates, N = 9,423
complete_pids <- spouse_df %>% 
  dplyr::filter(fup == 0) %>% 
  # delete rows with any of these missing
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(c(
        "dm", "fpg", "hba1c", "age", "sex", "edu_category", "employ_category", "hhincome",
        "bmi", "famhx_dm", "htn", "hld", "chd", "smk_overall", "alc_overall",
        "religion", "caste"
      )),
      ~ !is.na(.)
    )
  ) %>% 
  distinct(pid)

complete_df <- spouse_df %>% 
  dplyr::filter(pid %in% complete_pids$pid)


# Step 3: SPOUSE - ONE MALE + ONE FEMALE

# Define valid dyads: exactly 2 people in the same household, 1 male + 1 female 
valid_hhids <- complete_df %>%
  distinct(hhid, pid, sex) %>%       # one row per person
  group_by(hhid) %>%
  dplyr::filter(
    n() == 2,                        # exactly 2 people in the household
    all(sex %in% c("male","female")),           # only male/female
    sum(sex == "male") == 1,             # one male
    sum(sex == "female") == 1              # one female
  ) %>%
  ungroup() %>%
  pull(hhid) %>%
  unique()

# N = 7,196, 3,598 couples
analytic_spouses <- complete_df %>%
  dplyr::filter(hhid %in% valid_hhids) %>% 
  arrange(hhid, pid) 

# Step 4: Keep SPOUSE - AGE GAP <= 18Y 

# N = 3,518 spouses
age_gap18 <- analytic_spouses %>% 
  dplyr::filter(fup == 0) %>% 
  distinct(hhid, pid, age) %>%
  group_by(hhid) %>% 
  reframe(age_diff = abs(diff(age))) %>% 
  dplyr::filter(age_diff <= 18)

# exclude spouses with age gap >18y, N = 7,036, 3,518 couples
analytic_age18 <- analytic_spouses %>% 
  dplyr::filter(hhid %in% age_gap18$hhid) 

saveRDS(analytic_age18, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre01_long_spouse_complete_cases.RDS"))


# Archive: Derive incident diabetes events and time-to-event for baseline DM-free participants.
# Steps: Identify incident/censor dates and merge back to analytic data.
# Output: Saves psdcpre03_analytic dataset long.RDS.
rm(list=ls());gc();source(".Rprofile")

recoded_df <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre02_recoded_df.RDS"))

# Identify dm_biomarker free population at baseline, N = 2,029
dm_free_ids <- recoded_df %>% 
  dplyr::filter(fup == 0 & dm_biomarker == 0) %>% 
  distinct(pid)

# Identify dm_biomarker events at baseline, N = 41
dm_baseline_ids <- recoded_df %>% 
  dplyr::filter(fup == 0 & dm_biomarker != 0) %>% 
  distinct(pid)


# Incident dm_biomarker dates among baseline dm_biomarker-free participants only, N = 453
dm_incident_df <- recoded_df %>%
  semi_join(dm_free_ids, by = "pid") %>%
  arrange(carrs, pid, fup, doi) %>%
  group_by(carrs, pid) %>%
  dplyr::filter(dm_biomarker == 1) %>%
  slice_head(n = 1) %>%
  ungroup() %>% 
  mutate(dm_incident_date = inter_visit_midpoint,
         dm_incident = 1) %>% 
  mutate(time_to_dm = as.numeric(difftime(dm_incident_date, baseline_date, units = "days"))) %>% 
  select(carrs, fup, pid, hhid, dm_incident, dm_incident_date, time_to_dm) 

# Never have dm_biomarker, censored, N = 1,576
dm_never_df <- recoded_df %>%
  semi_join(dm_free_ids, by = "pid") %>%
  arrange(carrs, pid, fup, doi) %>%
  group_by(carrs, pid) %>%
  dplyr::filter(!any(dm_biomarker == 1, na.rm = TRUE)) %>%
  slice_tail(n = 1) %>%
  ungroup() %>% 
  mutate(dm_incident_date = doi,
         dm_incident = 0) %>% 
  mutate(time_to_dm = as.numeric(difftime(dm_incident_date, baseline_date, units = "days"))) %>% 
  select(carrs, fup, pid, hhid, dm_incident, dm_incident_date, time_to_dm) 

dm_event_or_censor_df <- dplyr::bind_rows(dm_incident_df, dm_never_df)
  
  
analytic_df <- recoded_df %>% 
  left_join(dm_event_or_censor_df,
            by = c("carrs", "fup", "pid", "hhid")) %>% 
  group_by(carrs, pid) %>% 
  mutate(baseline_dm = case_when(pid %in% dm_baseline_ids$pid ~ 1,
                                 TRUE ~ 0)) %>% 
  ungroup()

saveRDS(analytic_df, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre03_analytic dataset long.RDS"))


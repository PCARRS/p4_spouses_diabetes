# Archive: Convert long analytic dataset to spouse-wide format for dyad analyses.
# Steps: Pivot to female_/male_ columns and keep households with both spouses.
# Output: Saves psdcpre04_analytic dataset wide.RDS.
rm(list=ls());gc();source(".Rprofile")

analytic_df <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre03_analytic dataset long.RDS"))

# convert into “husband‑wife” wide format

value_cols <- setdiff(
  names(analytic_df),
  c("hhid", "sex", "carrs", "fup")    # drop your id‐cols here
)
# unique hhid: 1,035
analytic_df_wide <- analytic_df %>%
  pivot_wider(
    id_cols    = c(hhid, carrs, fup),   
    names_from = sex,
    values_from = all_of(value_cols),
    names_glue = "{sex}_{.value}"
  ) %>% 
  # Keep only household-wave records where both spouses are present after widening.
  dplyr::filter(!is.na(female_pid), !is.na(male_pid))


saveRDS(analytic_df_wide, paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/archive/psdcpre04_analytic_dataset_wide.RDS"))


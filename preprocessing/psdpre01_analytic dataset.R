rm(list=ls());gc();source(".Rprofile")

# Imputed wide dataset
harmonized_datasets <- readRDS(paste0(path_spouses_bmi_change_folder, "/working/cleaned/psban01_imputed harmonized dfs.RDS"))

spousedyads_clean <- readRDS(paste0(path_spouses_bmi_change_folder, "/working/preprocessing/spouseyads cleaned.RDS"))

wide_list <- vector("list", length(harmonized_datasets))

for (imp_num in seq_along(harmonized_datasets)) {
  df <- harmonized_datasets[[imp_num]]

  df1 <- df %>%
    mutate(visit = case_when(
      carrs == 1 & fup == 0 ~ 0,
      carrs == 1 & fup == 2 ~ 1,
      carrs == 1 & fup == 4 ~ 2,
      carrs == 1 & fup == 7 ~ 3,
      carrs == 2 & fup == 0 ~ 0,
      carrs == 2 & fup == 2 ~ 1,
      TRUE ~ NA_real_
    )) %>%
    mutate(
      visitcomplete = case_when(
        is.na(reason) ~ 1,
        TRUE ~ 0
      ),
      visit_date = case_when(
        visitcomplete == 1 ~ doi,
        TRUE ~ as.Date(NA)
      ),
      dm_biomarker = case_when(
        visitcomplete == 1 & is.na(fpg) & is.na(hba1c) ~ NA_real_,
        visitcomplete == 1 & (fpg >= 126 | hba1c >= 6.5 | dm == 1) ~ 1,
        visitcomplete == 1 ~ 0,
        TRUE ~ NA_real_
      ),
      dm_biomarker_nomiss = case_when(
        visitcomplete == 1 & (is.na(fpg) | is.na(hba1c)) ~ NA_real_,
        visitcomplete == 1 & (fpg >= 126 | hba1c >= 6.5 | dm == 1) ~ 1,
        visitcomplete == 1 ~ 0,
        TRUE ~ NA_real_
      ),
      dm_self = case_when(
        visitcomplete == 1 & dm == 1 ~ 1,
        visitcomplete == 1 ~ 0,
        TRUE ~ NA_real_
      ),
      missingbiomarkern = case_when(
        visitcomplete == 1 ~ rowSums(is.na(cbind(fpg, hba1c))),
        TRUE ~ NA_real_
      )
    )

  wide_outcomes <- df1 %>%
    select(
      pid, visit, fpg, hba1c, dm, dm_biomarker, dm_biomarker_nomiss,
      dm_self, visit_date, visitcomplete, missingbiomarkern
    ) %>%
    tidyr::pivot_wider(
      names_from = visit,
      values_from = c(
        fpg, hba1c, dm, dm_biomarker, dm_biomarker_nomiss,
        dm_self, visit_date, visitcomplete, missingbiomarkern
      ),
      names_glue = "{.value}{visit}"
    )

  baseline <- df1 %>%
    dplyr::filter(visit == 0) %>%
    dplyr::select(
      carrs, hhid, pid, doi, dob, ceb, age, sex, site, educstat, employ, hhincome,
      smk_overall, smk_exp, alc_overall, htn, dm, hld, chd, cva, ckd, cancer,
      famhx_htn, famhx_cvd, famhx_dm, sbp1, sbp2, sbp3, dbp1, dbp2, dbp3,
      height_cm, weight_kg, bmi, waist_cm, hip_cm
    )

  wide_list[[imp_num]] <- baseline %>%
    left_join(wide_outcomes, by = "pid") %>%
    left_join(
      spousedyads_clean %>% select(pid, hhid, spousedyad_new),
      by = c("pid", "hhid")
    )
}

saveRDS(
  wide_list,
  paste0(path_spouses_diabetes_folder, "/working/preprocessing/psdpre01_wide analytic dataset_list.RDS")
)



# Purpose: Descriptive characteristics for complete-case dyads.
# Steps: Compute sex-specific summaries, correlations, and missingness.
# Output: Writes psdcan01_descriptive characteristics.csv.
rm(list=ls());gc();source(".Rprofile")

dyads <- readRDS(paste0(path_spouses_diabetes_folder,"/working/cca/preprocessing/psdcpre03_spouse dyad dataset.RDS"))

continuous_vars <- c("age", "bmi", "waist_cm")  
proportion_vars <- c("smk_overall", "alc_overall", "famhx_htn", "famhx_dm", "famhx_cvd",
                     "chd", "cva", "ckd", "dm_biomarker0", "htn")
grouped_vars <- c("educcat", "employocccat",  "multimorbiditycat", "bmicat")     

# 1. Continuous variables: Mean (SD) and Pearson correlation
continuous_tbl <- map_dfr(continuous_vars, function(var) {
  female_col <- paste0(var, "_wife")
  male_col <- paste0(var, "_husb")
  
  f_vals <- dyads[[female_col]]
  m_vals <- dyads[[male_col]]
  
  # Calculate correlation for paired data
  r <- cor(f_vals, m_vals, use = "pairwise.complete.obs")
  
  tibble(
    variable = var,
    female = sprintf("%.1f (%.1f)", mean(f_vals, na.rm = TRUE), sd(f_vals, na.rm = TRUE)),
    male = sprintf("%.1f (%.1f)", mean(m_vals, na.rm = TRUE), sd(m_vals, na.rm = TRUE)),
    compare = sprintf("%.3f", r),
    missing_female = sprintf("%.1f%%", 100 * mean(is.na(f_vals))),
    missing_male = sprintf("%.1f%%", 100 * mean(is.na(m_vals)))
  )
})

# 2. Binary proportion variables
binary_tbl <- map_dfr(proportion_vars, function(var) {
  female_col <- paste0(var, "_wife")
  male_col <- paste0(var, "_husb")
  
  f_vals <- dyads[[female_col]]
  m_vals <- dyads[[male_col]]
  
  # Counts and percentages
  f_n <- sum(f_vals == 1, na.rm = TRUE)
  m_n <- sum(m_vals == 1, na.rm = TRUE)
  f_N <- sum(!is.na(f_vals))
  m_N <- sum(!is.na(m_vals))
  
  # Create 2x2 contingency table for OR calculation (male x female)
  tbl <- table(m_vals, f_vals)
  
  or_result <- tryCatch({
    if(nrow(tbl) == 2 && ncol(tbl) == 2) {
      fisher_result <- fisher.test(tbl)
      list(
        or = fisher_result$estimate,
        lower = fisher_result$conf.int[1],
        upper = fisher_result$conf.int[2]
      )
    } else {
      list(or = NA_real_, lower = NA_real_, upper = NA_real_)
    }
  }, error = function(e) {
    list(or = NA_real_, lower = NA_real_, upper = NA_real_)
  })
  
  # Format OR with CI
  if(is.na(or_result$or)) {
    or_ci <- "NA"
  } else {
    or_ci <- sprintf("%.2f (%.2f, %.2f)", or_result$or, or_result$lower, or_result$upper)
  }
  
  tibble(
    variable = var,
    female = sprintf("%d (%.1f%%)", f_n, 100 * f_n / f_N),
    male = sprintf("%d (%.1f%%)", m_n, 100 * m_n / m_N),
    compare = or_ci,
    missing_female = sprintf("%.1f%%", 100 * mean(is.na(f_vals))),
    missing_male = sprintf("%.1f%%", 100 * mean(is.na(m_vals)))
  )
})

# 3. Grouped variables
categorical_tbl <- map_dfr(grouped_vars, function(var) {
  female_col <- paste0(var, "_wife")
  male_col <- paste0(var, "_husb")
  
  f_vals <- dyads[[female_col]]
  m_vals <- dyads[[male_col]]
  
  # Get all levels
  all_levels <- unique(c(f_vals, m_vals))
  all_levels <- all_levels[!is.na(all_levels)]
  all_levels <- sort(all_levels)
  
  # Chi-square test
  chi_tab <- table(data.frame(
    var = c(f_vals, m_vals),
    sex = rep(c("female", "male"), times = c(length(f_vals), length(m_vals)))
  ))
  
  chi_result <- tryCatch({
    chisq.test(chi_tab)
  }, error = function(e) NULL)
  
  pval <- if(!is.null(chi_result)) chi_result$p.value else NA_real_
  pval_fmt <- if(is.na(pval)) "NA" else 
    if(pval < 0.001) "< 0.001" else 
      sprintf("%.3f", pval)
  
  f_miss <- 100 * mean(is.na(f_vals))
  m_miss <- 100 * mean(is.na(m_vals))
  
  # Calculate proportions for each level
  result_rows <- map_dfr(seq_along(all_levels), function(i) {
    level <- all_levels[i]
    f_prop <- sum(f_vals == level, na.rm = TRUE) / sum(!is.na(f_vals)) * 100
    m_prop <- sum(m_vals == level, na.rm = TRUE) / sum(!is.na(m_vals)) * 100
    
    tibble(
      variable = if(i == 1) var else "",
      level = as.character(level),
      female = sprintf("%.1f%%", f_prop),
      male = sprintf("%.1f%%", m_prop),
      compare = if(i == 1) pval_fmt else "",
      missing_female = if(i == 1) sprintf("%.1f%%", f_miss) else "",
      missing_male = if(i == 1) sprintf("%.1f%%", m_miss) else ""
    )
  })
  
  result_rows
})

# Add level column to continuous and binary tables
continuous_tbl <- continuous_tbl %>% mutate(level = NA_character_, .after = variable)
binary_tbl <- binary_tbl %>% mutate(level = NA_character_, .after = variable)

# Reorder columns consistently
continuous_tbl <- continuous_tbl %>% select(variable, level, female, male, compare, missing_female, missing_male)
binary_tbl <- binary_tbl %>% select(variable, level, female, male, compare, missing_female, missing_male)
categorical_tbl <- categorical_tbl %>% select(variable, level, female, male, compare, missing_female, missing_male)

# Combine all
table1 <- bind_rows(continuous_tbl, binary_tbl, categorical_tbl)


write.csv(table1, "cca/analysis/psdcan01_descriptive characteristics.csv", row.names = FALSE)

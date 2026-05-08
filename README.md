# Spousal Diabetes Study

Study objective: evaluate spousal concordance and longitudinal risk of diabetes in couples using baseline characteristics and follow-up outcomes.

The `cca/` folder contains complete-case analysis outputs; complete-case means only records with no missing values for the analysis variables are used. 
The top-level `analysis/` and `preprocessing/` folders use multiple imputation datasets rather than complete-case data.

## File guide

### Top-level

- `p4_spouses_diabetes.Rproj` is the RStudio project file for this workspace.
- `.Rprofile` stores R session options used when this project loads. Please add your folder path here.

### preprocessing/

- `preprocessing/psdpre01_analytic dataset.R` builds the analytic dataset from raw inputs.
- `preprocessing/psdpre02_recoded dataset.R` recodes variables and prepares analysis-ready fields.
- `preprocessing/psdpre03_spouse dyad dataset.R` constructs the spouse dyad dataset used in analyses.

### analysis/

- `analysis/psdan01_descriptive characteristics.R` generates descriptive characteristics summaries.
- `analysis/psdan02_dm freq logistic.R` runs diabetes concordance frequency and logistic models.
- `analysis/psdan03_followup risk by spouse baseline.R` computes follow-up risk by spouse baseline diabetes status.
- `analysis/psdan04_debug_single_imp.R` is a debugging script for single-imputation runs.
- `analysis/psdan04_spousal_baseline_dm_logistic.R` fits spousal baseline diabetes logistic models.
- `analysis/psdan05_spousal_baseline_dm_rr.R` fits spousal baseline diabetes risk ratio models.
- `analysis/psdan06_spouse_cox_models.R` fits spouse Cox models for incident diabetes.

### cca/preprocessing/

- `cca/preprocessing/psdcpre01_analytic dataset.R` builds the complete-case analytic dataset.
- `cca/preprocessing/psdcpre02_recoded dataset.R` recodes variables for complete-case analyses.
- `cca/preprocessing/psdcpre03_spouse dyad dataset.R` constructs complete-case spouse dyads.

### cca/analysis/

- `cca/analysis/psdcan01_descriptive characteristics.R` generates descriptive characteristics for complete cases.
- `cca/analysis/psdcan02_dm freq logistic.R` runs complete-case concordance frequency and logistic models.
- `cca/analysis/psdcan03_followup risk by spouse baseline.R` computes complete-case follow-up risk by spouse baseline diabetes status.
- `cca/analysis/psdcan04_spousal_baseline_dm_logistic.R` fits complete-case spousal baseline diabetes logistic models.
- `cca/analysis/psdcan05_spousal_baseline_dm_rr.R` fits complete-case spousal baseline diabetes risk ratio models.
- `cca/analysis/psdcan06_spouse_cox_models.R` fits complete-case spouse Cox models.

### cca/analysis/archive/ and cca/preprocessing/archive/

- Files in `cca/analysis/archive/` and `cca/preprocessing/archive/` are older or alternative versions kept for reference and are not part of the primary run order.

## How to run (order)

1. Run preprocessing scripts in order: `preprocessing/psdpre01_analytic dataset.R`, `preprocessing/psdpre02_recoded dataset.R`, `preprocessing/psdpre03_spouse dyad dataset.R`.
2. Run analysis scripts in order: `analysis/psdan01_descriptive characteristics.R` through `analysis/psdan06_spouse_cox_models.R`.
3. For complete-case analysis, repeat the same order under `cca/preprocessing/` then `cca/analysis/`.


Author:** Jiali Guo
**Email:** jguo2581@gmail.com
**Institution:** Emory University  
**Last Updated:** May 8, 2026

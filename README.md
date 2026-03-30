# Survival Modeling and Fairness Analysis in Metastatic NSCLC

This repository contains R scripts used to perform survival analysis, model evaluation, and fairness assessment in patients with metastatic non-small cell lung cancer (NSCLC) using SEER data.

The codebase is designed to ensure **transparency, reproducibility, and clarity** of the analytical workflow used in the study.

---

## Overview

The analysis includes:

- Cohort construction and preprocessing
- Descriptive statistics and baseline characteristics
- Kaplan–Meier survival analysis
- Cox proportional hazards modeling
- Machine learning survival models (RSF, BlackBoost)
- Model performance evaluation (C-index, IBS, iAUC)
- Calibration analysis
- Fairness assessment using the Keya et al. framework
- Intersectional performance analysis
- Visualization of results (forest plots, calibration curves, fairness plots)

---

## Repository Structure

All scripts are located in the root directory and are intended to be run sequentially.

### Data Preparation and Descriptive Analysis
- `01_data_preprocessing.R`  
  Cohort selection, data cleaning, and variable preparation.

- `02_descriptive_tables.R`  
  Generation of Table 1 (baseline characteristics) and Table 2 (alive/dead counts over time).

---

### Survival Analysis
- `03_km_subgroup_plots.R`  
  Kaplan–Meier survival curves across key subgroups.

- `04_cox_tables_full_and_train.R`  
  Univariable and multivariable Cox regression analyses (Table 3 and Table 4).

---

### Data Splitting
- `05_data_split.R`  
  Stratified train–test split based on event status.

---

### Model Evaluation
- `06_main_model_evaluation.R`  
  Comprehensive model evaluation pipeline including performance metrics and alternative fairness calculations (used for exploratory analysis).

- `07_metrics_train_test.R`  
  Primary performance evaluation (C-index, IBS, iAUC) on training and test datasets.

---

### Fairness Analysis
- `08_fairness_keya.R`  
  Main fairness analysis based on the framework proposed by Keya et al., including:
  - Individual fairness (Fi)
  - Group fairness (Fg)
  - Concordance imparity (CI%)
  - Intersectional fairness

---

### Advanced Analyses
- `09_intersectional_performance.R`  
  Intersectional subgroup performance analysis to identify heterogeneity in model discrimination.

- `10_model_comparison_and_timing.R`  
  Supplementary analysis including model comparison and computational benchmarking.

---

### Calibration Analysis
- `11_calibration_subgroups.R`  
  Subgroup-specific calibration curves at multiple time horizons (24, 36, 60 months).

- `12_calibration_12m_all_models.R`  
  Calibration analysis at a fixed 12-month horizon across all models.

---

### Visualization Scripts
- `13_forestplot_multivariable_cox.R`  
  Forest plot for multivariable Cox model with hazard ratios and event counts.

- `14_fairness_barplot.R`  
  Bar plots summarizing fairness metrics (CI%, Fi, Fg) across models and subgroups.

- `15_lollipop_subgroup_performance.R`  
  Lollipop plots illustrating subgroup-specific performance (C-index, IBS, iAUC).

- `16_patient_selection_flowchart.R`  
  Flowchart describing cohort selection and exclusion criteria (Figure 1).

---

## Reproducibility

- All analyses are conducted in **R**
- Random processes use fixed seeds for reproducibility
- Scripts are modular and can be executed independently after preprocessing

---

## Notes

- The primary results reported in the manuscript are based on:
  - `07_metrics_train_test.R`
  - `08_fairness_keya.R`
  - `09_intersectional_performance.R`

- Additional scripts (e.g., model benchmarking and alternative evaluation pipelines) are provided for completeness and transparency.

---

## Disclaimer

This repository is intended for research purposes and reflects the analytical workflow used in the associated study.

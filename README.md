# Survival Modeling and Fairness Analysis in Metastatic NSCLC

This repository contains R scripts used to perform survival analysis, model evaluation, and fairness assessment in patients with metastatic non-small cell lung cancer (NSCLC) using SEER data.

The codebase is designed to ensure **transparency, reproducibility, and clarity** of the analytical workflow used in the study.

---

## Table of Contents

1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Environment Setup](#environment-setup)
4. [Input Data Format](#input-data-format)
5. [Models](#models)
6. [Performance Evaluation](#performance-evaluation)
7. [How to Run the Code](#how-to-run-the-code)
8. [Data Availability](#data-availability)
9. [Reproducibility](#reproducibility)
10. [Notes](#notes)
11. [Disclaimer](#disclaimer)

---

## Overview

The analysis includes:

- Cohort construction and preprocessing  
- Descriptive statistics and baseline characteristics  
- Kaplan-Meier survival analysis  
- Cox proportional hazards modeling  
- Machine learning survival models (RSF, BlackBoost)  
- Model performance evaluation (C-index, IBS, iAUC)  
- Calibration analysis  
- Fairness assessment using the Keya et al. framework  
- Intersectional performance analysis  
- Visualization of results (forest plots, calibration curves, fairness plots)

---

## Repository Structure

All scripts are located in the repository root and are intended to be run sequentially.

### Data Preparation and Descriptive Analysis

- `01_data_preprocessing.R`  
  Performs cohort selection, applies inclusion and exclusion criteria, and prepares the final analytic dataset.

- `02_descriptive_tables.R`  
  Generates:
  - Table 1: baseline characteristics  
  - Table 2: alive/dead counts over fixed time points  

### Survival Analysis

- `03_km_subgroup_plots.R`  
  Produces Kaplan–Meier survival curves stratified by key demographic and clinical subgroups.

- `04_cox_tables_full_and_train.R`  
  Performs univariable and multivariable Cox regression analyses and generates:
  - Table 3: full dataset results  
  - Table 4: training dataset results  

### Data Splitting

- `05_data_split.R`  
  Performs a stratified 70/30 train–test split based on event status.

### Model Evaluation

- `06_main_model_evaluation.R`  
  Comprehensive evaluation pipeline including model fitting, performance estimation, and exploratory fairness calculations.  
  This script is retained for methodological completeness.

- `07_metrics_train_test.R`  
  Primary performance evaluation used in the manuscript.  
  Computes:
  - C-index  
  - Integrated Brier Score (IBS)  
  - Integrated AUC (iAUC)  

### Fairness Analysis

- `08_fairness_keya.R`  
  Implements fairness metrics based on the framework proposed by Keya et al., including:
  - Individual fairness (Fi)  
  - Group fairness (Fg)  
  - Concordance imparity (CI%)  
  - Intersectional fairness  

### Advanced Analyses

- `09_intersectional_performance.R`  
  Evaluates model performance across intersectional subgroup combinations.

- `10_model_comparison_and_timing.R`  
  Provides supplementary analyses including model comparison and computational benchmarking.

### Calibration Analysis

- `11_calibration_subgroups.R`  
  Generates subgroup-specific calibration curves at multiple time horizons.

- `12_calibration_12m_all_models.R`  
  Generates calibration plots at a fixed 12-month horizon.

### Visualization Scripts

- `13_forestplot_multivariable_cox.R`  
  Produces the forest plot for the multivariable Cox model.

- `14_fairness_barplot.R`  
  Generates bar plots summarizing fairness metrics.

- `15_lollipop_subgroup_performance.R`  
  Produces lollipop plots illustrating subgroup-specific performance.

- `16_patient_selection_flowchart.R`  
  Generates the cohort selection flowchart.

---

## Environment Setup

All analyses were conducted in **R (version 4.5.0)**.

### Required R Packages

- survival  
- survminer  
- randomForestSRC  
- mboost  
- pec  
- riskRegression  
- dplyr  
- tidyverse  
- ggplot2  
- patchwork  
- cowplot  
- DiagrammeR  
- DiagrammeRsvg  
- rsvg  
- caret  
- broom  
- stringr  

Install packages using:

    install.packages(c(
      "survival","survminer","randomForestSRC","mboost",
      "pec","riskRegression","dplyr","tidyverse",
      "ggplot2","patchwork","cowplot","DiagrammeR",
      "DiagrammeRsvg","rsvg","caret","broom","stringr"
    ))

---

## Input Data Format

The analysis uses a dataset derived from the SEER database.

### Raw Data Source

The analytic dataset is derived from:

- `seer_nsclc_metastatic_dataset.txt`

This file represents the original SEER extract used as input for preprocessing.

### Main Input File

- `data_n19254.csv` — final analytic cohort (n = 19,254)

### Data Processing

The preprocessing pipeline implemented in `01_data_preprocessing.R`:

- applies inclusion and exclusion criteria  
- constructs derived variables  
- reproduces the cohort definition described in Wang et al.  

### Required Variables

**Time-to-event variables**
- Survival.months  
- Vital.status.recode  

**Demographic variables**
- Age_grouped  
- Gender  
- Race_grouped  
- Marital_grouped  
- Income_group  
- Urban_Rural  

**Clinical variables**
- Primary_site_grouped  
- Histology_grouped  
- Grade_grouped  
- Laterality_grouped  
- T_stage_recoded  
- N_stage_combined  
- M_stage_combined  

**Treatment variables**
- Surgery  
- Chemotherapy  
- Radiation  

**Metastasis indicators**
- Bone_met  
- Liver_met  
- Lung_met  
- Brain_met  

### Notes

- `event_status` is generated within scripts:
  - 1 = death  
  - 0 = censored  

---

## Models

- Cox Proportional Hazards (CoxPH)  
- Random Survival Forest (RSF)  
- BlackBoost  

---

## Performance Evaluation

Model performance and fairness were evaluated using complementary metrics.

### Discrimination and Prediction Error

- C-index  
- Integrated Brier Score (IBS)  
- Integrated time-dependent AUC (iAUC)  

### Fairness Metrics

- Individual Fairness (Fi)  
- Group Fairness (Fg)  
- Intersectional Fairness (F∩)  
- Concordance Imparity (CI%)  

### Evaluation Strategy

- Metrics were computed on both training and test datasets  
- Performance and fairness were assessed across demographic and clinical subgroups  

---

## How to Run the Code

Run scripts sequentially:

1. `01_data_preprocessing.R`  
2. `02_descriptive_tables.R`  
3. `03_km_subgroup_plots.R`  
4. `04_cox_tables_full_and_train.R`  
5. `05_data_split.R`  
6. `07_metrics_train_test.R`  
7. `08_fairness_keya.R`  
8. `09_intersectional_performance.R`  

---

## Data Availability

The data used in this study were obtained from the SEER database.

Due to data use agreements, the dataset cannot be publicly shared.

https://seer.cancer.gov/

---

## Reproducibility

- Fixed random seeds are used  
- Scripts are modular  

---

## Notes

Primary results are based on:

- `07_metrics_train_test.R`  
- `08_fairness_keya.R`  
- `09_intersectional_performance.R`  

---

## Disclaimer

This repository is intended for research purposes only.

# Survival Fairness Metrics

This repository contains R scripts developed for survival analysis and comparative evaluation of model performance and fairness.  
The codebase supports methodological analyses across multiple survival modeling approaches and subgroup structures.

The repository is intended to provide a transparent and reproducible record of analytical workflows used in ongoing research.

## Repository Structure

All analysis scripts are located in the repository root and are designed to be executed sequentially.

- **1.Data_Selection.R** – Cohort selection and outcome definition  
- **2.Table1&2.R** – Descriptive statistics and baseline tables  
- **3.KM_Curves.R** – Kaplan–Meier survival curves  
- **4.Cox_Uni_Multi.R** – Univariable and multivariable Cox regression models  
- **5.DataSplit.R** – Train–test data splitting procedure  
- **6.Comparisons.R** – Model-level and subgroup comparisons  
- **6.Metrics_TrainTest.R** – Performance metrics on training and test sets  
- **7.CC.R** – Calibration curve analysis  
- **8.CT.R** – Comparison table generation  
- **9.Int_Pairwise.R** – Intersectional pairwise analyses  
- **CC-12.R** – Fixed-horizon (12-month) analysis  
- **Flowchart.R** – Cohort flowchart generation  
- **Keya.R** – Fairness metrics based on published methodology  
- **barchart.R** – Bar chart visualizations  
- **lollipop.R** – Lollipop plot visualizations  


# =========================================================
# Script: 02_descriptive_tables.R
# Description:
# Generates Table 1 (descriptive summary) and Table 2
# (alive/dead counts over time) from the cleaned dataset.
#
# Input:
# - data_n19254.csv
#
# Output:
# - summary_table_n19254.csv
# - table2_alive_dead_over_time.csv
# =========================================================

# =========================
# Load packages
# =========================
library(dplyr)
library(readr)

# =========================
# Load data
# =========================
data <- read.csv("data_n19254.csv")

# =========================
# Create event indicator
# =========================
# 1 = dead, 0 = alive
data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)

# =========================
# TABLE 1: Descriptive summary
# =========================
vars_to_summarize <- c(
  "Age_grouped",
  "Race_grouped",
  "Gender",
  "Marital_grouped",
  "Urban_Rural",
  "Income_group",
  "Primary_site_grouped",
  "Histology_grouped",
  "Grade_grouped",
  "Laterality_grouped",
  "T_stage_recoded",
  "N_stage_combined",
  "M_stage_combined",
  "Surgery",
  "Chemotherapy",
  "Radiation",
  "Bone_met",
  "Liver_met",
  "Lung_met",
  "Brain_met",
  "event_status"
)

summarize_variable <- function(varname) {
  if (!varname %in% names(data)) {
    warning(paste("Variable not found:", varname))
    return(NULL)
  }
  
  tab <- table(data[[varname]], useNA = "no")
  
  if (length(tab) == 0) {
    warning(paste("No data for variable:", varname))
    return(NULL)
  }
  
  df <- data.frame(
    Variable = varname,
    Category = names(tab),
    Count = as.vector(tab),
    Percentage = round(as.vector(tab) / sum(tab) * 100, 1)
  )
  
  return(df)
}

summary_table <- do.call(rbind, lapply(vars_to_summarize, summarize_variable))

# Save Table 1
write_excel_csv(summary_table, "summary_table_n19254.csv")

# =========================
# TABLE 2: Alive/dead counts over time
# =========================
# Cumulative deaths are calculated at fixed time points
# using the full cohort as the denominator.

total_n <- nrow(data)

time_points <- c(12, 24, 36, 48, 60)

status_fixed_total <- lapply(time_points, function(tp) {
  dead <- sum(data$Survival.months <= tp & data$event_status == 1, na.rm = TRUE)
  alive <- total_n - dead
  
  data.frame(
    Time_Month = tp,
    Alive = alive,
    Alive_Pct = round(100 * alive / total_n, 1),
    Dead = dead,
    Dead_Pct = round(100 * dead / total_n, 1),
    Total = total_n
  )
})

time_status_fixed_total <- do.call(rbind, status_fixed_total)

# Save Table 2
write_excel_csv(time_status_fixed_total, "table2_alive_dead_over_time.csv")

cat("Table 1 and Table 2 outputs saved.\n")
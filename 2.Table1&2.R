
### TABLE 1 ###

data <- read.csv("data_n19254.csv")
# Define variables to summarize
vars_to_summarize <- c(
  "Age_grouped",
  "Race_grouped",
  "Gender",
  "Marital_grouped",
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
  "Urban_Rural",
  "Income_group"
)

# Safe summary function: handles missing or absent variables gracefully
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

# Apply the function to all selected variables and combine results
summary_table <- do.call(rbind, lapply(vars_to_summarize, summarize_variable))

# Save the summary table as CSV
write.csv(summary_table, "summary_table_n19254.csv", row.names = FALSE, fileEncoding = "UTF-8")


library(readr)

write_excel_csv(summary_table, "summary_table_n19254.csv")

# Load packages
library(dplyr)
library(readr)

# Read data
data <- read.csv("data_n19254.csv")
# Create event indicator (1 = dead, 0 = alive)
data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)
# ==== TABLE 1 ====
vars_to_summarize <- c(
  "Age_grouped", "Race_grouped", "Gender", "Marital_grouped", "Urban_Rural", "Income_group",
  "Primary_site_grouped","Histology_grouped", "Grade_grouped", "Laterality_grouped", "T_stage_recoded",
  "N_stage_combined", "M_stage_combined", "Surgery", "Chemotherapy", "Radiation",
  "Bone_met", "Liver_met", "Lung_met", "Brain_met", "event_status"
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

# ==== Alive/Dead Counts at Time Points ====
# ---- 3. Cumulative Deaths Over Time (Fixed Total) ----
library(dplyr)
library(readr)

# Define total number of patients
total_n <- nrow(data)  # should be 19254

# Define time points
time_points <- c(12, 24, 36, 48, 60)

# Calculate cumulative deaths at each time point
status_fixed_total <- lapply(time_points, function(tp) {
  dead <- sum(data$Survival.months <= tp & data$event_status == 1, na.rm = TRUE)
  alive <- total_n - dead
  
  alive_pct <- round(100 * alive / total_n, 1)
  dead_pct <- round(100 * dead / total_n, 1)
  
  data.frame(
    Time_Month = tp,
    Alive = alive,
    Alive_Pct = alive_pct,
    Dead = dead,
    Dead_Pct = dead_pct,
    Total = total_n
  )
})

# Combine and save
time_status_fixed_total <- do.call(rbind, status_fixed_total)
write_excel_csv(time_status_fixed_total, "alive_dead_fixed_total.csv")
cat("✅ Time-based cumulative deaths with fixed total saved.\n")

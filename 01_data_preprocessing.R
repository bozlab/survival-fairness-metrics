# =========================================================
# Project: Survival and Fairness Analysis in Metastatic NSCLC
# Description:
# This script prepares the analytic dataset from SEER data
# for replication of Wang et al. and subsequent modeling.
#
# Input:
# - seer_nsclc_metastatic_dataset.txt
#
# Output:
# - data_n19254.csv
# - variable_summary_table.csv
# =========================================================

# =========================
# Load required libraries
# =========================
library(dplyr)

# =========================
# Load dataset
# =========================
# SEER-derived tab-delimited dataset for metastatic NSCLC
data <- read.delim("seer_nsclc_metastatic_dataset.txt",
                   na.strings = "NA",
                   check.names = TRUE)

cat("Initial number of patients diagnosed with lung cancer:", nrow(data), "\n")

# =========================
# Cohort selection
# =========================

# 1. Filter by diagnosis year
data <- subset(data, Year.of.diagnosis >= 2004 & Year.of.diagnosis <= 2016)
cat("After filtering by diagnosis year (2004–2016):", nrow(data), "\n")

# 2. Define histology groups
squamous_cell_codes <- c(8051:8052, 8070:8076, 8078, 8083:8084, 8090, 8094, 8123)
adenocarcinoma_codes <- c(
  8015, 8050, 8140:8141, 8143:8145, 8147, 8190, 8201, 8211,
  8250:8255, 8260, 8290, 8310, 8320, 8323, 8333, 8401, 8440,
  8470:8471, 8480:8481, 8490, 8503, 8507, 8550, 8570:8572,
  8574, 8576, 8012:8014, 8021, 8034, 8082
)
large_cell_codes <- c(
  8046, 8003:8004, 8022, 8030:8033, 8035, 8120, 8200,
  8240:8241, 8243:8246, 8249, 8430, 8525, 8560, 8562, 8575
)
small_cell_codes <- c(8002, 8041:8045)

# 3. Add histology group
data$Histology_Group <- with(
  data,
  ifelse(Histologic.Type.ICD.O.3 %in% squamous_cell_codes, "Squamous Cell Carcinoma",
         ifelse(Histologic.Type.ICD.O.3 %in% adenocarcinoma_codes, "Adenocarcinoma",
                ifelse(Histologic.Type.ICD.O.3 %in% large_cell_codes, "Large Cell Carcinoma",
                       ifelse(Histologic.Type.ICD.O.3 %in% small_cell_codes, "Small Cell Lung Cancer",
                              "Other/Unspecified"))))
)
cat("After defining histology groups:", nrow(data), "\n")

# 4. Exclude small cell lung cancer
data <- subset(data, Histology_Group != "Small Cell Lung Cancer")
cat("After excluding small cell histology (NSCLC only):", nrow(data), "\n")

# 5. Keep only M1 stage (metastatic)
data <- subset(data, grepl("^M1", Derived.AJCC.M..7th.ed..2010.2015.))
cat("After selecting M1 stage only (exclude M0, MX, unknown):", nrow(data), "\n")

# 6. Exclude patients with unknown metastatic site data
data <- subset(
  data,
  SEER.Combined.Mets.at.DX.bone..2010.. %in% c("Yes", "No") &
    SEER.Combined.Mets.at.DX.brain..2010.. %in% c("Yes", "No") &
    SEER.Combined.Mets.at.DX.liver..2010.. %in% c("Yes", "No") &
    SEER.Combined.Mets.at.DX.lung..2010.. %in% c("Yes", "No")
)
cat("After excluding unknown metastatic sites (bone, brain, liver, lung):", nrow(data), "\n")

# 7. Exclude patients under 18 years
data <- subset(data, !grepl("^0|1-4|5-9|10-14|15-19", Age.recode.with..1.year.olds))
cat("After excluding patients under 18 years:", nrow(data), "\n")

# 8. Keep only first malignant primary tumors
data <- subset(data, First.malignant.primary.indicator == "Yes")
cat("After filtering for first malignant primary tumors:", nrow(data), "\n")

# 9. Exclude patients with >1 malignant tumor
data <- subset(data, Total.number.of.in.situ.malignant.tumors.for.patient == "01")
cat("After excluding patients with >1 malignant tumor:", nrow(data), "\n")

# 10. Exclude patients with missing, unknown, or zero survival months
data$Survival.months <- suppressWarnings(as.numeric(as.character(data$Survival.months)))
data <- subset(data, !is.na(Survival.months) & Survival.months > 0)
cat("After removing cases with missing, unknown, or 0 months survival time:", nrow(data), "\n")

# =========================
# Covariate cleaning
# =========================

# 11. Recode Laterality into grouped categories
data$Laterality_grouped <- dplyr::case_when(
  data$Laterality == "Left - origin of primary" ~ "Left",
  data$Laterality == "Right - origin of primary" ~ "Right",
  data$Laterality == "Bilateral, single primary" ~ "Bilateral",
  TRUE ~ NA_character_
)

# 12. Exclude patients with missing/unknown covariates and surgery info
data <- subset(
  data,
  !is.na(Race.recode..White..Black..Other.) &
    Race.recode..White..Black..Other. != "Unknown" &
    !is.na(Marital.status.at.diagnosis) &
    Marital.status.at.diagnosis != "Unknown" &
    !is.na(Primary.Site...labeled) &
    Primary.Site...labeled != "C34.9-Lung, NOS" &
    !is.na(Grade.Recode..thru.2017.) &
    Grade.Recode..thru.2017. != "Unknown" &
    !is.na(Laterality_grouped) &
    !is.na(RX.Summ..Surg.Prim.Site..1998..) &
    RX.Summ..Surg.Prim.Site..1998.. != 99
)
cat("After removing missing/unknown covariates (and C34.9):", nrow(data), "\n")

# 13. Exclude TX and NX stage cases
data <- subset(
  data,
  !(
    (!is.na(Derived.AJCC.T..7th.ed..2010.2015.) &
       Derived.AJCC.T..7th.ed..2010.2015. == "TX") |
      (!is.na(Derived.AJCC.T..6th.ed..2004.2015.) &
         Derived.AJCC.T..6th.ed..2004.2015. == "TX")
  ) &
    !(
      (!is.na(Derived.AJCC.N..7th.ed..2010.2015.) &
         Derived.AJCC.N..7th.ed..2010.2015. == "NX") |
        (!is.na(Derived.AJCC.N..6th.ed..2004.2015.) &
           Derived.AJCC.N..6th.ed..2004.2015. == "NX")
    )
)
cat("After excluding TX and NX stage cases:", nrow(data), "\n")

# =========================
# Summary table
# =========================

vars_to_summarize <- c(
  "Age.recode.with..1.year.olds",
  "Sex",
  "Race.recode..White..Black..Other.",
  "Marital.status.at.diagnosis",
  "Primary.Site...labeled",
  "Histology_Group",
  "Grade.Recode..thru.2017.",
  "Laterality",
  "RX.Summ..Surg.Prim.Site..1998..",
  "Radiation.recode",
  "Chemotherapy.recode..yes..no.unk.",
  "SEER.Combined.Mets.at.DX.bone..2010..",
  "SEER.Combined.Mets.at.DX.brain..2010..",
  "SEER.Combined.Mets.at.DX.liver..2010..",
  "SEER.Combined.Mets.at.DX.lung..2010..",
  "Derived.AJCC.T..7th.ed..2010.2015.",
  "Derived.AJCC.N..7th.ed..2010.2015.",
  "Derived.AJCC.M..7th.ed..2010.2015.",
  "Derived.AJCC.T..6th.ed..2004.2015.",
  "Derived.AJCC.N..6th.ed..2004.2015.",
  "Derived.AJCC.M..6th.ed..2004.2015.",
  "Derived.AJCC.Stage.Group..6th.ed..2004.2015.",
  "Vital.status.recode..study.cutoff.used."
)

summarize_variable <- function(varname) {
  tab <- table(data[[varname]])
  total <- sum(tab)
  df <- data.frame(
    Variable = varname,
    Category = names(tab),
    Count = as.integer(tab),
    Percent = round(100 * as.integer(tab) / total, 1)
  )
  df$Summary <- paste0(df$Count, " (", df$Percent, "%)")
  df[, c("Variable", "Category", "Summary")]
}

summary_table <- do.call(rbind, lapply(vars_to_summarize, summarize_variable))
write.csv(summary_table, "variable_summary_table.csv", row.names = FALSE)
cat("Summary table saved as 'variable_summary_table.csv'.\n")

# =========================
# Stage combination and recoding
# =========================

combine_T_stage <- function(row) {
  if (!is.na(row["Derived.AJCC.T..7th.ed..2010.2015."]) &&
      row["Derived.AJCC.T..7th.ed..2010.2015."] != "TX") {
    return(row["Derived.AJCC.T..7th.ed..2010.2015."])
  } else if (!is.na(row["Derived.AJCC.T..6th.ed..2004.2015."]) &&
             row["Derived.AJCC.T..6th.ed..2004.2015."] != "TX") {
    return(row["Derived.AJCC.T..6th.ed..2004.2015."])
  } else {
    return(NA)
  }
}
data$T_stage_combined <- apply(data, 1, combine_T_stage)

recode_T <- function(stage) {
  if (stage %in% c("T0", "T1", "T1a", "T1b", "T1NOS")) return("T0+T1")
  if (stage %in% c("T2", "T2a", "T2b", "T2NOS")) return("T2")
  if (stage == "T3") return("T3")
  if (stage == "T4") return("T4")
  return(NA)
}
data$T_stage_recoded <- sapply(data$T_stage_combined, recode_T)

combine_N_stage <- function(row) {
  if (!is.na(row["Derived.AJCC.N..7th.ed..2010.2015."]) &&
      row["Derived.AJCC.N..7th.ed..2010.2015."] != "NX") {
    return(row["Derived.AJCC.N..7th.ed..2010.2015."])
  } else if (!is.na(row["Derived.AJCC.N..6th.ed..2004.2015."]) &&
             row["Derived.AJCC.N..6th.ed..2004.2015."] != "NX") {
    return(row["Derived.AJCC.N..6th.ed..2004.2015."])
  } else {
    return(NA)
  }
}
data$N_stage_combined <- apply(data, 1, combine_N_stage)

combine_M_stage <- function(row) {
  if (!is.na(row["Derived.AJCC.M..7th.ed..2010.2015."])) {
    return(row["Derived.AJCC.M..7th.ed..2010.2015."])
  } else if (!is.na(row["Derived.AJCC.M..6th.ed..2004.2015."])) {
    return(row["Derived.AJCC.M..6th.ed..2004.2015."])
  } else {
    return(NA)
  }
}
data$M_stage_combined <- apply(data, 1, combine_M_stage)

# =========================
# Analysis variable recoding
# =========================

data$Surgery <- ifelse(data$RX.Summ..Surg.Prim.Site..1998.. == "0", "No", "Yes")

data$Chemotherapy <- ifelse(
  data$Chemotherapy.recode..yes..no.unk. == "Yes",
  "Yes",
  "No/unknown"
)

data$Radiation <- dplyr::case_when(
  data$Radiation.recode %in% c(
    "Beam radiation",
    "Combination of beam with implants or isotopes",
    "Radioactive implants (includes brachytherapy) (1988+)",
    "Radioisotopes (1988+)",
    "Radiation, NOS  method or source not specified"
  ) ~ "Yes",
  TRUE ~ "No"
)

data$Bone_met  <- data$SEER.Combined.Mets.at.DX.bone..2010..
data$Liver_met <- data$SEER.Combined.Mets.at.DX.liver..2010..
data$Lung_met  <- data$SEER.Combined.Mets.at.DX.lung..2010..
data$Brain_met <- data$SEER.Combined.Mets.at.DX.brain..2010..

data$Race_grouped <- ifelse(
  data$Race.recode..White..Black..Other. == "White", "White",
  ifelse(data$Race.recode..White..Black..Other. == "Black", "Black", "Other")
)

data$Marital_grouped <- ifelse(
  data$Marital.status.at.diagnosis == "Married (including common law)",
  "Married",
  "Others"
)

data$Age_grouped <- case_when(
  data$Age.recode.with..1.year.olds %in% c(
    "20-24 years", "25-29 years", "30-34 years", "35-39 years",
    "40-44 years", "45-49 years", "50-54 years"
  ) ~ "20–54 years",
  data$Age.recode.with..1.year.olds %in% c("55-59 years", "60-64 years") ~ "55–64 years",
  data$Age.recode.with..1.year.olds %in% c("65-69 years", "70-74 years") ~ "65–74 years",
  data$Age.recode.with..1.year.olds %in% c("75-79 years", "80-84 years") ~ "75–84 years",
  data$Age.recode.with..1.year.olds == "85+ years" ~ "85+ years",
  TRUE ~ NA_character_
)

data$Gender <- data$Sex

data$Histology_grouped <- ifelse(
  data$Histology_Group == "Adenocarcinoma", "Adenocarcinoma",
  ifelse(data$Histology_Group == "Squamous Cell Carcinoma", "Squamous", "Others")
)

data$Grade_grouped <- ifelse(
  data$Grade.Recode..thru.2017. == "Well differentiated; Grade I", "Grade I",
  ifelse(data$Grade.Recode..thru.2017. == "Moderately differentiated; Grade II", "Grade II",
         ifelse(data$Grade.Recode..thru.2017. == "Poorly differentiated; Grade III", "Grade III",
                ifelse(data$Grade.Recode..thru.2017. == "Undifferentiated; anaplastic; Grade IV",
                       "Grade IV", NA)))
)

data$Primary_site_grouped <- ifelse(
  data$Primary.Site...labeled == "C34.0-Main bronchus", "Main bronchus",
  ifelse(data$Primary.Site...labeled == "C34.1-Upper lobe, lung", "Upper lobe",
         ifelse(data$Primary.Site...labeled == "C34.2-Middle lobe, lung", "Middle lobe",
                ifelse(data$Primary.Site...labeled == "C34.3-Lower lobe, lung", "Lower lobe",
                       ifelse(data$Primary.Site...labeled == "C34.8-Overlapping lesion of lung",
                              "Overlapping lesion of lung", NA))))
)

data$Urban_Rural <- ifelse(
  data$Rural.Urban.Continuum.Code %in% c(
    "Counties in metropolitan areas ge 1 million pop",
    "Counties in metropolitan areas of 250,000 to 1 million pop",
    "Counties in metropolitan areas of lt 250 thousand pop"
  ),
  "Urban",
  ifelse(
    data$Rural.Urban.Continuum.Code %in% c(
      "Nonmetropolitan counties adjacent to a metropolitan area",
      "Nonmetropolitan counties not adjacent to a metropolitan area"
    ),
    "Rural",
    NA
  )
)

data$Income_group <- case_when(
  data$Median.household.income.inflation.adj.to.2022 %in% c(
    "< $40,000", "$40,000 - $44,999", "$45,000 - $49,999", "$50,000 - $54,999"
  ) ~ "< $55,000",
  data$Median.household.income.inflation.adj.to.2022 %in% c(
    "$55,000 - $59,999", "$60,000 - $64,999", "$65,000 - $69,999", "$70,000 - $74,999"
  ) ~ "$55,000 - $75,000",
  data$Median.household.income.inflation.adj.to.2022 %in% c(
    "$75,000 - $79,999", "$80,000 - $84,999", "$85,000 - $89,999"
  ) ~ "$75,000 - $89,999",
  data$Median.household.income.inflation.adj.to.2022 %in% c(
    "$90,000 - $94,999", "$95,000 - $99,999", "$100,000 - $109,999",
    "$110,000 - $119,999", "$120,000+"
  ) ~ "> $90,000",
  TRUE ~ NA_character_
)

# =========================
# Save final cleaned dataset
# =========================
write.csv(data, "data_n19254.csv", row.names = FALSE)
cat("Final cleaned dataset saved as 'data_n19254.csv'.\n")

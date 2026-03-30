# Load required libraries
library(survival)
library(dplyr)
library(broom)
library(readr)

# Load full dataset
data <- read.csv("data_n19254.csv")

# Create event_status if not already present
if (!"event_status" %in% names(data)) {
  data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)
}

# Define survival object
surv_obj <- Surv(time = data$Survival.months, event = data$event_status)

# Define predictor variables
predictors <- c("Age_grouped", "Race_grouped", "Gender", "Marital_grouped",
                "Primary_site_grouped", "Histology_grouped", "Grade_grouped",
                "Laterality_grouped", "T_stage_recoded", "N_stage_combined",
                "M_stage_combined", "Surgery", "Chemotherapy", "Radiation",
                "Bone_met", "Liver_met", "Lung_met", "Brain_met")

# ---------- 1. Univariable Cox regression ----------
uni_results <- lapply(predictors, function(var) {
  formula <- as.formula(paste("surv_obj ~", var))
  model <- coxph(formula, data = data)
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(variable = var)
})
uni_df <- bind_rows(uni_results) %>%
  rename(HR = estimate, CI_lower = conf.low, CI_upper = conf.high, p = p.value) %>%
  select(variable, term, HR, CI_lower, CI_upper, p)

# ---------- 2. Multivariable Cox regression ----------
significant_vars <- unique(uni_df$variable[uni_df$p < 0.05])
multi_formula <- as.formula(paste("surv_obj ~", paste(significant_vars, collapse = " + ")))
multi_model <- coxph(multi_formula, data = data)
multi_summary <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(HR = estimate, CI_lower = conf.low, CI_upper = conf.high, p = p.value) %>%
  select(term, HR, CI_lower, CI_upper, p)

# ---------- 3. Event summary function ----------
get_event_summary <- function(varname, dataset) {
  if (!varname %in% names(dataset)) return(NULL)
  df <- dataset %>%
    group_by(!!sym(varname)) %>%
    summarise(
      Event = sum(event_status == 1, na.rm = TRUE),
      Total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Event_Pct = paste0(Event, " (", round(100 * Event / Total, 1), "%)"),
      term = paste0(varname, !!sym(varname))
    ) %>%
    select(term, Event_Pct)
  return(df)
}
event_summary <- bind_rows(lapply(predictors, get_event_summary, dataset = data))

# ---------- 4. Format HR (CI) and p-values ----------
format_hr_ci <- function(HR, lower, upper) {
  sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)
}
format_p <- function(pvec) {
  sapply(pvec, function(pval) {
    if (is.na(pval)) return("")
    if (pval < 0.001) "<.001" else sprintf("%.3f", pval)
  })
}
uni_df_fmt <- uni_df %>%
  mutate(
    HR_CI_uni = format_hr_ci(HR, CI_lower, CI_upper),
    p_uni = format_p(p)
  ) %>%
  select(term, variable, HR_CI_uni, p_uni)

multi_df_fmt <- multi_summary %>%
  mutate(
    HR_CI_multi = format_hr_ci(HR, CI_lower, CI_upper),
    p_multi = format_p(p)
  ) %>%
  select(term, HR_CI_multi, p_multi)

# ---------- 5. Merge all ----------
table3 <- full_join(uni_df_fmt, multi_df_fmt, by = "term") %>%
  left_join(event_summary, by = "term")

# ---------- 6. Add reference rows WITH Event_Pct ----------
get_reference_row <- function(varname, dataset) {
  ref_level <- levels(as.factor(dataset[[varname]]))[1]
  ref_data <- dataset %>% filter(!!sym(varname) == ref_level)
  
  # Event and total count for reference level
  event_n <- sum(ref_data$event_status == 1, na.rm = TRUE)
  total_n <- nrow(ref_data)
  event_pct <- paste0(event_n, " (", round(100 * event_n / total_n, 1), "%)")
  
  term <- paste0(varname, ref_level)
  data.frame(
    term = term,
    variable = varname,
    HR_CI_uni = "Reference",
    p_uni = "–",
    HR_CI_multi = "Reference",
    p_multi = "–",
    Event_Pct = event_pct,
    stringsAsFactors = FALSE
  )
}

ref_rows <- bind_rows(lapply(predictors, get_reference_row, dataset = data))

# Remove duplicates before binding (in case already present)
table3 <- bind_rows(table3, ref_rows) %>%
  distinct(term, .keep_all = TRUE)
# ---------- 7. Sort with reference at top within each group ----------
table3_final <- bind_rows(lapply(predictors, function(var) {
  df_var <- table3 %>% filter(variable == var)
  df_var <- df_var %>%
    mutate(is_ref = ifelse(HR_CI_uni == "Reference", 1, 0)) %>%
    arrange(desc(is_ref), term) %>%
    select(-is_ref)
  return(df_var)
}))

# ---------- 8. Final selection and save ----------
table3_final <- table3_final %>%
  select(term, Event_Pct, HR_CI_uni, p_uni, HR_CI_multi, p_multi)

write.csv(table3_final, "Table3_final.csv", row.names = FALSE)
cat("✅ Final Table 3 saved with reference categories inline.\n")




##### Train Data #####

# Load required libraries
library(survival)
library(dplyr)
library(broom)
library(readr)

# Load full dataset
data <- read.csv("train_data_seer.csv")

# Create event_status if not already present
if (!"event_status" %in% names(data)) {
  data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)
}

# Define survival object
surv_obj <- Surv(time = data$Survival.months, event = data$event_status)

# Define predictor variables
predictors <- c("Age_grouped", "Race_grouped", "Gender", "Marital_grouped",
                "Primary_site_grouped", "Histology_grouped", "Grade_grouped",
                "Laterality_grouped", "T_stage_recoded", "N_stage_combined",
                "M_stage_combined", "Surgery", "Chemotherapy", "Radiation",
                "Bone_met", "Liver_met", "Lung_met", "Brain_met")

# ---------- 1. Univariable Cox regression ----------
uni_results <- lapply(predictors, function(var) {
  formula <- as.formula(paste("surv_obj ~", var))
  model <- coxph(formula, data = data)
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(variable = var)
})
uni_df <- bind_rows(uni_results) %>%
  rename(HR = estimate, CI_lower = conf.low, CI_upper = conf.high, p = p.value) %>%
  select(variable, term, HR, CI_lower, CI_upper, p)

# ---------- 2. Multivariable Cox regression ----------
significant_vars <- unique(uni_df$variable[uni_df$p < 0.05])
multi_formula <- as.formula(paste("surv_obj ~", paste(significant_vars, collapse = " + ")))
multi_model <- coxph(multi_formula, data = data)
multi_summary <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(HR = estimate, CI_lower = conf.low, CI_upper = conf.high, p = p.value) %>%
  select(term, HR, CI_lower, CI_upper, p)

# ---------- 3. Event summary function ----------
get_event_summary <- function(varname, dataset) {
  if (!varname %in% names(dataset)) return(NULL)
  df <- dataset %>%
    group_by(!!sym(varname)) %>%
    summarise(
      Event = sum(event_status == 1, na.rm = TRUE),
      Total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Event_Pct = paste0(Event, " (", round(100 * Event / Total, 1), "%)"),
      term = paste0(varname, !!sym(varname))
    ) %>%
    select(term, Event_Pct)
  return(df)
}
event_summary <- bind_rows(lapply(predictors, get_event_summary, dataset = data))

# ---------- 4. Format HR (CI) and p-values ----------
format_hr_ci <- function(HR, lower, upper) {
  sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)
}
format_p <- function(pvec) {
  sapply(pvec, function(pval) {
    if (is.na(pval)) return("")
    if (pval < 0.001) "<.001" else sprintf("%.3f", pval)
  })
}
uni_df_fmt <- uni_df %>%
  mutate(
    HR_CI_uni = format_hr_ci(HR, CI_lower, CI_upper),
    p_uni = format_p(p)
  ) %>%
  select(term, variable, HR_CI_uni, p_uni)

multi_df_fmt <- multi_summary %>%
  mutate(
    HR_CI_multi = format_hr_ci(HR, CI_lower, CI_upper),
    p_multi = format_p(p)
  ) %>%
  select(term, HR_CI_multi, p_multi)

# ---------- 5. Merge all ----------
table3 <- full_join(uni_df_fmt, multi_df_fmt, by = "term") %>%
  left_join(event_summary, by = "term")

# ---------- 6. Add reference rows WITH Event_Pct ----------
get_reference_row <- function(varname, dataset) {
  ref_level <- levels(as.factor(dataset[[varname]]))[1]
  ref_data <- dataset %>% filter(!!sym(varname) == ref_level)
  
  # Event and total count for reference level
  event_n <- sum(ref_data$event_status == 1, na.rm = TRUE)
  total_n <- nrow(ref_data)
  event_pct <- paste0(event_n, " (", round(100 * event_n / total_n, 1), "%)")
  
  term <- paste0(varname, ref_level)
  data.frame(
    term = term,
    variable = varname,
    HR_CI_uni = "Reference",
    p_uni = "–",
    HR_CI_multi = "Reference",
    p_multi = "–",
    Event_Pct = event_pct,
    stringsAsFactors = FALSE
  )
}

ref_rows <- bind_rows(lapply(predictors, get_reference_row, dataset = data))

# Remove duplicates before binding (in case already present)
table3 <- bind_rows(table3, ref_rows) %>%
  distinct(term, .keep_all = TRUE)
# ---------- 7. Sort with reference at top within each group ----------
table3_final <- bind_rows(lapply(predictors, function(var) {
  df_var <- table3 %>% filter(variable == var)
  df_var <- df_var %>%
    mutate(is_ref = ifelse(HR_CI_uni == "Reference", 1, 0)) %>%
    arrange(desc(is_ref), term) %>%
    select(-is_ref)
  return(df_var)
}))

# ---------- 8. Final selection and save ----------
table3_final <- table3_final %>%
  select(term, Event_Pct, HR_CI_uni, p_uni, HR_CI_multi, p_multi)

write.csv(table3_final, "Table4_train.csv", row.names = FALSE)




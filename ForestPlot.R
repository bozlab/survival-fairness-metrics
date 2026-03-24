# ============================================================
# FULL FINAL FOREST PLOT
# Multivariable Cox model with:
# Variable | Events / N (%) | Forest | Adjusted HR (95% CI)
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(broom)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# 1. Load data
# ============================================================
data <- read.csv("train_data_seer.csv", stringsAsFactors = FALSE)

if (!"event_status" %in% names(data)) {
  data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)
}

data <- data %>%
  filter(!is.na(Survival.months), !is.na(event_status))

# ============================================================
# 2. Predictors
# ============================================================
predictors <- c(
  "Age_grouped", "Race_grouped", "Gender", "Marital_grouped",
  "Primary_site_grouped", "Histology_grouped", "Grade_grouped",
  "Laterality_grouped", "T_stage_recoded", "N_stage_combined",
  "M_stage_combined", "Surgery", "Chemotherapy", "Radiation",
  "Bone_met", "Liver_met", "Lung_met", "Brain_met"
)

predictors <- predictors[predictors %in% names(data)]
data[predictors] <- lapply(data[predictors], as.factor)

# ============================================================
# 2b. Manual factor ordering for selected variables
# ============================================================
if ("Race_grouped" %in% names(data)) {
  current_levels <- levels(as.factor(data$Race_grouped))
  wanted_levels  <- c("Black", "White", "Other")
  wanted_levels  <- wanted_levels[wanted_levels %in% current_levels]
  other_levels   <- setdiff(current_levels, wanted_levels)
  data$Race_grouped <- factor(
    data$Race_grouped,
    levels = c(wanted_levels, other_levels)
  )
  data$Race_grouped <- droplevels(data$Race_grouped)
}

if ("Histology_grouped" %in% names(data)) {
  current_levels <- levels(as.factor(data$Histology_grouped))
  wanted_levels  <- c("Adenocarcinoma", "Squamous", "Others")
  wanted_levels  <- wanted_levels[wanted_levels %in% current_levels]
  other_levels   <- setdiff(current_levels, wanted_levels)
  data$Histology_grouped <- factor(
    data$Histology_grouped,
    levels = c(wanted_levels, other_levels)
  )
  data$Histology_grouped <- droplevels(data$Histology_grouped)
}

# ============================================================
# 3. Survival object
# ============================================================
surv_obj <- Surv(time = data$Survival.months, event = data$event_status)

# ============================================================
# 4. Univariable screening
# ============================================================
uni_results <- lapply(predictors, function(var) {
  formula_uni <- as.formula(paste("surv_obj ~", var))
  model_uni <- coxph(formula_uni, data = data)
  
  tidy(model_uni, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(variable = var)
})

uni_df <- bind_rows(uni_results)
significant_vars <- unique(uni_df$variable[uni_df$p.value < 0.05])

if (length(significant_vars) == 0) {
  stop("No significant predictors were identified in univariable screening.")
}

# ============================================================
# 5. Force variable order in multivariable model
# ============================================================
desired_order <- c(
  "Age_grouped",
  "Gender",
  "Race_grouped",
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
  "Brain_met"
)

significant_vars <- desired_order[desired_order %in% significant_vars]

# ============================================================
# 6. Multivariable model
# ============================================================
multi_formula <- as.formula(
  paste("surv_obj ~", paste(significant_vars, collapse = " + "))
)

multi_model <- coxph(multi_formula, data = data)

# ============================================================
# 7. Variable label map
# ============================================================
var_label_map <- c(
  Age_grouped          = "Age",
  Race_grouped         = "Race",
  Gender               = "Sex",
  Marital_grouped      = "Marital status",
  Primary_site_grouped = "Primary site",
  Histology_grouped    = "Histology",
  Grade_grouped        = "Grade",
  Laterality_grouped   = "Laterality",
  T_stage_recoded      = "T stage",
  N_stage_combined     = "N stage",
  M_stage_combined     = "M stage",
  Surgery              = "Surgery",
  Chemotherapy         = "Chemotherapy",
  Radiation            = "Radiation",
  Bone_met             = "Bone metastasis",
  Liver_met            = "Liver metastasis",
  Lung_met             = "Lung metastasis",
  Brain_met            = "Brain metastasis"
)

get_var_label <- function(var) {
  if (var %in% names(var_label_map)) {
    unname(var_label_map[var])
  } else {
    var
  }
}

# ============================================================
# 8. Tidy multivariable model
# ============================================================
tidy_multi <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(
    HR = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high,
    p = p.value
  )

# ============================================================
# 9. Build grouped forest dataframe from model$xlevels
# ============================================================
make_forest_df <- function(model, tidy_df, desired_order) {
  xlv <- model$xlevels
  vars_in_model <- desired_order[desired_order %in% names(xlv)]
  out <- list()
  
  for (var in vars_in_model) {
    levs <- xlv[[var]]
    ref_lev <- levs[1]
    nonref_levs <- levs[-1]
    var_lab <- get_var_label(var)
    
    header <- data.frame(
      variable = var,
      variable_label = var_lab,
      level = var_lab,
      HR = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      p = NA_real_,
      is_reference = NA,
      row_type = "header",
      row_order = 0,
      stringsAsFactors = FALSE
    )
    
    ref_row <- data.frame(
      variable = var,
      variable_label = var_lab,
      level = ref_lev,
      HR = 1,
      CI_lower = 1,
      CI_upper = 1,
      p = NA_real_,
      is_reference = TRUE,
      row_type = "body",
      row_order = 1,
      stringsAsFactors = FALSE
    )
    
    nonref_rows <- lapply(seq_along(nonref_levs), function(i) {
      lev <- nonref_levs[i]
      
      idx <- which(
        str_starts(tidy_df$term, fixed(var)) &
          str_ends(tidy_df$term, fixed(lev))
      )
      
      if (length(idx) == 0) return(NULL)
      
      rr <- tidy_df[idx[1], , drop = FALSE]
      
      data.frame(
        variable = var,
        variable_label = var_lab,
        level = lev,
        HR = rr$HR,
        CI_lower = rr$CI_lower,
        CI_upper = rr$CI_upper,
        p = rr$p,
        is_reference = FALSE,
        row_type = "body",
        row_order = i + 1,
        stringsAsFactors = FALSE
      )
    })
    
    out[[var]] <- bind_rows(header, ref_row, bind_rows(nonref_rows))
  }
  
  bind_rows(out)
}

plot_df <- make_forest_df(multi_model, tidy_multi, desired_order)

# ============================================================
# 10. Event counts per level
# ============================================================
get_event_counts <- function(data, var) {
  data %>%
    group_by(.data[[var]]) %>%
    summarise(
      events = sum(event_status == 1, na.rm = TRUE),
      N = n(),
      pct = 100 * events / N,
      .groups = "drop"
    ) %>%
    rename(level = .data[[var]])
}

event_list <- lapply(names(multi_model$xlevels), function(v) {
  tmp <- get_event_counts(data, v)
  tmp$variable <- v
  tmp
})

event_df <- bind_rows(event_list)

plot_df <- plot_df %>%
  left_join(event_df, by = c("variable", "level"))

# ============================================================
# 11. Safe variable ordering in plot_df
# ============================================================
plot_df <- plot_df %>%
  mutate(
    group_order = match(variable, desired_order)
  ) %>%
  arrange(group_order, row_order)

# ============================================================
# 12. Pretty text
# ============================================================
plot_df <- plot_df %>%
  mutate(
    display_label = ifelse(row_type == "header", level, paste0("   ", level)),
    display_label = str_replace_all(display_label, "_", " "),
    event_text = ifelse(
      row_type == "header",
      "",
      sprintf("%d / %d (%.1f%%)", events, N, pct)
    ),
    hrci_text = ifelse(
      row_type == "header",
      "",
      ifelse(
        is_reference,
        "Reference",
        sprintf("%.2f (%.2f–%.2f)", HR, CI_lower, CI_upper)
      )
    )
  )

# Top-to-bottom order: first row should appear at the TOP
plot_df$y <- seq(nrow(plot_df), 1, by = -1)
y_top <- max(plot_df$y) + 1.4

# ============================================================
# 13. Axis limits for forest panel
# ============================================================
xmin <- min(plot_df$CI_lower, na.rm = TRUE)
xmax <- max(plot_df$CI_upper, na.rm = TRUE)

xmin <- min(0.25, floor(xmin * 10) / 10)
xmax <- max(4, ceiling(xmax * 10) / 10)

forest_breaks <- c(0.25, 0.5, 1, 2, 4, 8)
forest_breaks <- forest_breaks[forest_breaks >= xmin & forest_breaks <= xmax]

# ============================================================
# 14. Left labels panel
# ============================================================
p_labels <- ggplot(plot_df, aes(y = y, x = 0, label = display_label)) +
  geom_text(
    data = subset(plot_df, row_type == "header"),
    hjust = 0, fontface = "bold", size = 5.0
  ) +
  geom_text(
    data = subset(plot_df, row_type == "body"),
    hjust = 0, size = 4.6
  ) +
  annotate(
    "text",
    x = 0, y = y_top,
    label = "Variable",
    hjust = 0, fontface = "bold", size = 5.2
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.5, y_top)) +
  theme_void() +
  theme(
    plot.margin = margin(8, 0, 8, 10)
  )

# ============================================================
# 15. Event column panel
# ============================================================
p_event <- ggplot(plot_df, aes(y = y, x = 0, label = event_text)) +
  geom_text(hjust = 0, size = 4.6) +
  annotate(
    "text",
    x = 0, y = y_top,
    label = "Events / N (%)",
    hjust = 0, fontface = "bold", size = 5.2
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.5, y_top)) +
  theme_void() +
  theme(
    plot.margin = margin(8, 0, 8, 0)
  )

# ============================================================
# 16. Forest panel
# ============================================================
p_forest <- ggplot(subset(plot_df, row_type == "body"), aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_segment(
    data = subset(plot_df, row_type == "body" & !is_reference),
    aes(x = CI_lower, xend = CI_upper, yend = y),
    linewidth = 0.85
  ) +
  geom_point(
    data = subset(plot_df, row_type == "body" & !is_reference),
    aes(x = HR),
    size = 3.2
  ) +
  geom_point(
    data = subset(plot_df, row_type == "body" & is_reference),
    aes(x = HR),
    size = 3.2, shape = 1, stroke = 1
  ) +
  scale_x_log10(
    limits = c(xmin, xmax),
    breaks = forest_breaks
  ) +
  scale_y_continuous(limits = c(0.5, y_top)) +
  labs(x = "Adjusted hazard ratio", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 2, 8, 2),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold")
  )

# ============================================================
# 17. HR (CI) text panel
# ============================================================
p_hr <- ggplot(plot_df, aes(y = y, x = 0, label = hrci_text)) +
  geom_text(hjust = 0, size = 4.6) +
  annotate(
    "text",
    x = 0, y = y_top,
    label = "Adjusted HR (95% CI)",
    hjust = 0, fontface = "bold", size = 5.2
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.5, y_top)) +
  theme_void() +
  theme(
    plot.margin = margin(8, 10, 8, 0)
  )

# ============================================================
# 18. Combine panels
# ============================================================
final_plot <- p_labels + p_event + p_forest + p_hr +
  plot_layout(widths = c(2.2, 1.8, 3.3, 2.3))

# ============================================================
# 19. Save
# ============================================================
plot_height <- max(8.5, 0.31 * nrow(plot_df))

ggsave(
  filename = "ForestPlot_Multivariable_Cox.png",
  plot = final_plot,
  width = 15,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "ForestPlot_Multivariable_Cox.pdf",
  plot = final_plot,
  width = 15,
  height = plot_height,
  bg = "white"
)

cat("✅ Final forest plot saved successfully.\n")

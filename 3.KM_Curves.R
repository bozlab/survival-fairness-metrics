# ============================================================
# Libraries
# ============================================================
library(survival)
library(survminer)
library(gridExtra)
library(dplyr)
library(forcats)   # fct_recode, fct_relevel

# ============================================================
# Data
# ============================================================
data <- read.csv("data_n19254.csv")

# Event indicator
data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)

# Keep valid rows
data <- subset(data, !is.na(Survival.months) & !is.na(event_status))

# Grouping variables
group_vars <- c("Age_grouped", "Gender", "Race_grouped", 
                "Marital_grouped", "Income_group", "Urban_Rural")

# Ensure factors
data[group_vars] <- lapply(data[group_vars], factor)

# ============================================================
# Factor level ordering / label fixes
# ============================================================

## ---- Race: Black, White, Other ----
data$Race_grouped <- factor(
  data$Race_grouped,
  levels = c("Black", "White", "Other")
)

## ---- Income: sadece LEVEL adını düzelt + sıralama ----
# "$55,000–$75,000" -> "$55,000–$74,999"
data$Income_group <- fct_recode(
  data$Income_group,
  "$55,000–$74,999" = "$55,000–$75,000"
)

# Küçükten büyüğe sıra: <55k, 55–74,999, 75–89,999, >90k
data$Income_group <- fct_relevel(
  data$Income_group,
  "< $55,000",
  "$55,000–$74,999",
  "$75,000–$89,999",
  "> $90,000"
)

# ============================================================
# Plot titles / legends
# ============================================================
plot_titles <- c(
  "Age",
  "Sex",
  "Race",
  "Marital Status",
  "Income",
  "Urban/Rural"
)

legend_titles <- c("Age", "Sex", "Race", "Marital Status", "Income", "Urban/Rural")

# Legend labels: only used levels, in specified order
legend_labels <- lapply(
  group_vars,
  function(var) levels(droplevels(data[[var]]))
)
# ---- Income legend manual order ----
legend_labels[[5]] <- c(
  "< $55,000",
  "$55,000–$74,999",
  "$75,000–$89,999",
  "> $90,000"
)

# ============================================================
# KM models
# ============================================================
km_models <- list(
  Age         = survfit(Surv(Survival.months, event_status) ~ Age_grouped,     data = data),
  Gender      = survfit(Surv(Survival.months, event_status) ~ Gender,          data = data),
  Race        = survfit(Surv(Survival.months, event_status) ~ Race_grouped,    data = data),
  Marital     = survfit(Surv(Survival.months, event_status) ~ Marital_grouped, data = data),
  Income      = survfit(Surv(Survival.months, event_status) ~ Income_group,    data = data),
  Urban_Rural = survfit(Surv(Survival.months, event_status) ~ Urban_Rural,     data = data)
)

# ============================================================
# KM plots (no risk table)
# ============================================================
km_plots <- mapply(
  function(model, title, legend_title, legend_lab) {
    p <- do.call(ggsurvplot, list(
      fit          = model,
      data         = data,
      title        = title,
      legend.title = legend_title,
      legend.labs  = legend_lab,
      risk.table   = FALSE,   # no risk table
      pval         = TRUE,
      pval.size    = 3,
      conf.int     = TRUE,
      xlab         = "Time in months",
      ylab         = "Survival probability",
      xlim         = c(0, 60),
      break.time.by = 12,
      palette      = "Dark2",
      surv.median.line = "hv",
      ggtheme = theme_minimal() +
        theme(
          plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title   = element_text(size = 8),
          axis.text    = element_text(size = 7),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8, face = "bold"),
          axis.line    = element_line(color = "black", linewidth = 0.8),
          plot.margin  = margin(5, 40, 5, 40)
        )
    ))
    
    # expand = 0 fix
    if (!"x" %in% names(p$plot$scales$scales)) {
      p$plot <- p$plot + scale_x_continuous(expand = c(0, 0))
    }
    if (!"y" %in% names(p$plot$scales$scales)) {
      p$plot <- p$plot + scale_y_continuous(expand = c(0, 0))
    }
    
    p$plot  # only the main plot (no risk table)
  },
  km_models,
  plot_titles,
  legend_titles,
  legend_labels,
  SIMPLIFY = FALSE
)

# ============================================================
# Save 2x3 grid: PNG + vector PDF
# ============================================================
png("KM_6plots_final_norisktable.png",
    width = 10, height = 10, units = "in", res = 300)
grid.arrange(grobs = km_plots, ncol = 2, nrow = 3)
dev.off()

pdf("~KM_6plots_final_norisktable.pdf",
    width = 10, height = 10)
grid.arrange(grobs = km_plots, ncol = 2, nrow = 3)
dev.off()

cat("KM plots (no risk table) saved successfully (PNG + PDF).\n")

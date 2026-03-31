# ============================================================
# Packages
# ============================================================
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# ============================================================
# 1) Subgroup performance results for all three models
#    (point estimates only; CIs stay in the Supplement table)
#    + Overall TEST SET row added at the top
# ============================================================
perf_raw <- tribble(
  ~Subgroup,    ~Level,              ~Model,       ~Cindex, ~IBS,   ~iAUC,
  
  #----------------- Overall (TEST SET) -----------------
  "Overall",    "Test set",          "CoxPH",        0.687,   0.065,  0.802,
  "Overall",    "Test set",          "RSF",          0.653,   0.063,  0.819,
  "Overall",    "Test set",          "BlackBoost",   0.694,   0.064,  0.815,
  
  #----------------- Sex -----------------
  "Sex",        "Female",            "CoxPH",        0.695,   0.082,  0.786,
  "Sex",        "Male",              "CoxPH",        0.680,   0.054,  0.809,
  
  "Sex",        "Female",            "RSF",          0.665,   0.079,  0.815,
  "Sex",        "Male",              "RSF",          0.643,   0.053,  0.822,
  
  "Sex",        "Female",            "BlackBoost",   0.705,   0.081,  0.809,
  "Sex",        "Male",              "BlackBoost",   0.686,   0.053,  0.819,
  
  #----------------- Age -----------------
  "Age",        "20–54",             "CoxPH",        0.671,   0.090,  0.769,
  "Age",        "55–64",             "CoxPH",        0.686,   0.070,  0.805,
  "Age",        "65–74",             "CoxPH",        0.696,   0.063,  0.783,
  "Age",        "75–84",             "CoxPH",        0.677,   0.052,  0.818,
  "Age",        "85+",               "CoxPH",        0.647,   0.064,  0.769,
  
  "Age",        "20–54",             "RSF",          0.643,   0.087,  0.798,
  "Age",        "55–64",             "RSF",          0.648,   0.067,  0.805,
  "Age",        "65–74",             "RSF",          0.657,   0.062,  0.792,
  "Age",        "75–84",             "RSF",          0.636,   0.051,  0.854,
  "Age",        "85+",               "RSF",          0.615,   0.061,  0.843,
  
  "Age",        "20–54",             "BlackBoost",   0.686,   0.088,  0.796,
  "Age",        "55–64",             "BlackBoost",   0.688,   0.068,  0.808,
  "Age",        "65–74",             "BlackBoost",   0.702,   0.062,  0.801,
  "Age",        "75–84",             "BlackBoost",   0.683,   0.052,  0.827,
  "Age",        "85+",               "BlackBoost",   0.649,   0.063,  0.808,
  
  #----------------- Race -----------------
  "Race",       "White",             "CoxPH",        0.687,   0.063,  0.805,
  "Race",       "Black",             "CoxPH",        0.695,   0.057,  0.702,
  "Race",       "Other",             "CoxPH",        0.688,   0.104,  0.742,
  
  "Race",       "White",             "RSF",          0.650,   0.062,  0.820,
  "Race",       "Black",             "RSF",          0.665,   0.054,  0.805,
  "Race",       "Other",             "RSF",          0.672,   0.099,  0.735,
  
  "Race",       "White",             "BlackBoost",   0.692,   0.063,  0.818,
  "Race",       "Black",             "BlackBoost",   0.700,   0.054,  0.746,
  "Race",       "Other",             "BlackBoost",   0.703,   0.100,  0.742,
  
  #----------------- Marital -----------------
  "Marital Status",    "Married",    "CoxPH",        0.690,   0.071,  0.799,
  "Marital Status",    "Others",     "CoxPH",        0.680,   0.058,  0.803,
  
  "Marital Status",    "Married",    "RSF",          0.653,   0.070,  0.809,
  "Marital Status",    "Others",     "RSF",          0.647,   0.055,  0.828,
  
  "Marital Status",    "Married",    "BlackBoost",   0.696,   0.071,  0.810,
  "Marital Status",    "Others",     "BlackBoost",   0.688,   0.057,  0.817,
  
  #----------------- Income -----------------
  "Income",     "< $55,000",         "CoxPH",        0.689,   0.065,  0.800,
  "Income",     "$55,000–$74,999",   "CoxPH",        0.687,   0.070,  0.781,
  "Income",     "$75,000–$89,999",   "CoxPH",        0.700,   0.069,  0.815,
  "Income",     "> $90,000",         "CoxPH",        0.680,   0.083,  0.781,
  
  "Income",     "< $55,000",         "RSF",          0.644,   0.065,  0.798,
  "Income",     "$55,000–$74,999",   "RSF",          0.658,   0.067,  0.807,
  "Income",     "$75,000–$89,999",   "RSF",          0.670,   0.067,  0.825,
  "Income",     "> $90,000",         "RSF",          0.649,   0.079,  0.793,
  
  "Income",     "< $55,000",         "BlackBoost",   0.693,   0.065,  0.796,
  "Income",     "$55,000–$74,999",   "BlackBoost",   0.696,   0.068,  0.801,
  "Income",     "$75,000–$89,999",   "BlackBoost",   0.704,   0.068,  0.827,
  "Income",     "> $90,000",         "BlackBoost",   0.685,   0.081,  0.791,
  
  #----------------- Urban/Rural -----------------
  "Urban/Rural","Urban",             "CoxPH",        0.691,   0.068,  0.806,
  "Urban/Rural","Rural",             "CoxPH",        0.678,   0.065,  0.783,
  
  "Urban/Rural","Urban",             "RSF",          0.658,   0.066,  0.818,
  "Urban/Rural","Rural",             "RSF",          0.641,   0.063,  0.809,
  
  "Urban/Rural","Urban",             "BlackBoost",   0.698,   0.067,  0.819,
  "Urban/Rural","Rural",             "BlackBoost",   0.683,   0.064,  0.793
)

# ---- GLOBAL LEVEL ORDER ----
perf_raw <- perf_raw %>%
  mutate(
    Level = factor(
      Level,
      levels = c(
        # Overall
        "Test set",
        # Sex
        "Female", "Male",
        # Age 
        "20–54", "55–64", "65–74", "75–84", "85+",
        # Race
        "Black", "White", "Other",
        # Marital
        "Married", "Others",
        # Income 
        "< $55,000", "$55,000–$74,999",
        "$75,000–$89,999", "> $90,000",
        # Urban / Rural 
        "Rural", "Urban"
      )
    )
  )

# ============================================================
# 2) Prepare long-format data and reference lines
# ============================================================
perf_long <- perf_raw %>%
  mutate(
    Subgroup = recode(Subgroup, "Urban / Rural" = "Urban/Rural"),
    Subgroup = factor(
      Subgroup,
      levels = c("Overall", "Sex", "Age", "Race",
                 "Marital Status", "Income", "Urban/Rural")
    ),
    Model = factor(Model, levels = c("CoxPH", "RSF", "BlackBoost"))
  ) %>%
  pivot_longer(
    cols = c(Cindex, IBS, iAUC),
    names_to = "Metric",
    values_to = "Value"
  )

line_data <- perf_long %>%
  group_by(Subgroup, Metric, Level) %>%
  summarise(
    xmin = case_when(
      Metric == "Cindex" ~ 0.50,
      Metric == "IBS"    ~ 0.00,
      Metric == "iAUC"   ~ 0.65
    ),
    xmax = max(Value),
    .groups = "drop"
  )

model_cols <- c(
  "CoxPH"      = "#4E79A7",
  "RSF"        = "#F28E2B",
  "BlackBoost" = "#59A14F"
)

# ============================================================
# 3) Helper function for one metric panel (no legend inside)
# ============================================================
make_lolli <- function(metric_name,
                       x_limits,
                       x_breaks,
                       x_lab,
                       show_y_labels = TRUE,
                       show_strips  = FALSE,
                       show_legend  = FALSE) {
  
  df_m   <- perf_long %>% filter(Metric == metric_name)
  line_m <- line_data %>% filter(Metric == metric_name)
  
  p <- ggplot() +
    geom_segment(
      data = line_m,
      aes(x = xmin, xend = xmax, y = Level, yend = Level),
      color = "grey85", linewidth = 0.4
    ) +
    geom_point(
      data = df_m,
      aes(x = Value, y = Level, color = Model),
      size = 2
    ) +
    facet_grid(
      Subgroup ~ .,
      scales = "free_y",
      space  = "free_y",
      labeller = labeller(
        Subgroup = c(
          "Overall"        = "Overall",
          "Sex"            = "Sex",
          "Age"            = "Age",
          "Race"           = "Race",
          "Marital Status" = "Marital\nStatus",
          "Income"         = "Income",
          "Urban/Rural"    = "Urban/Rural"
        )
      )
    ) +
    scale_color_manual(values = model_cols, name = "Model") +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    scale_y_discrete(limits = function(x) rev(x)) +
    labs(x = x_lab, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 11),
      axis.text.y        = element_text(size = 11),
      legend.position    = if (show_legend) "top" else "none",
      legend.direction   = "horizontal",
      legend.title       = element_text(face = "bold", size = 11),
      legend.text        = element_text(size = 10),
      plot.margin        = margin(t = 5, r = 5, b = 5, l = 5),
      strip.clip         = "off"
    )
  
  if (!show_y_labels) {
    p <- p +
      theme(
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  if (!show_strips) {
    p <- p +
      theme(
        strip.background = element_blank(),
        strip.text.y     = element_blank()
      )
  } else {
    p <- p +
      theme(
        strip.text.y = element_text(
          face   = "bold",
          angle  = 0,
          margin = margin(t = 2, b = 2, l = 2, r = 12),
          size   = 11
        )
      )
  }
  
  p
}

# ============================================================
# 4) Build the three metric panels
# ============================================================
p_cindex <- make_lolli(
  metric_name    = "Cindex",
  x_limits       = c(0.50, 0.72),
  x_breaks       = seq(0.50, 0.75, by = 0.05),
  x_lab          = "C-index",
  show_y_labels  = TRUE,
  show_strips    = FALSE,
  show_legend    = TRUE
)

p_ibs <- make_lolli(
  metric_name    = "IBS",
  x_limits       = c(0.00, 0.11),
  x_breaks       = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10),
  x_lab          = "IBS",
  show_y_labels  = FALSE,
  show_strips    = FALSE,
  show_legend    = FALSE
)

p_iauc <- make_lolli(
  metric_name    = "iAUC",
  x_limits       = c(0.65, 0.87),
  x_breaks       = seq(0.65, 0.85, by = 0.05),
  x_lab          = "iAUC",
  show_y_labels  = FALSE,
  show_strips    = TRUE,
  show_legend    = FALSE
)

# ============================================================
# 5) Combine panels & legend
# ============================================================
final_plot <- (p_cindex + p_ibs + p_iauc) +
  plot_layout(nrow = 1, guides = "collect") &
  theme(
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.title     = element_text(face = "bold")
  )

final_plot

# ============================================================
# 6) Save as PNG and vector PDF
# ============================================================
ggsave(
  "Figure_subgroup_performance_lollipop_3models_final_with_overall.png",
  final_plot, width = 11, height = 5, dpi = 300
)

ggsave(
  "Figure_subgroup_performance_lollipop_3models_final_with_overall.pdf",
  final_plot, width = 11, height = 5, device = "pdf"
)

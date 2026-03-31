library(tidyverse)

fair_raw <- tribble(
  ~Subgroup,     ~Model,       ~CI,    ~CI_low, ~CI_high, ~Fi,    ~Fi_low, ~Fi_high, ~Fg,    ~Fg_low, ~Fg_high,
  "Age",         "Cox",        5.980,  3.649,   8.903,    2.194,  2.154,   2.245,    1.127,  0.956,   1.527,
  "Age",         "RSF",        5.946,  3.038,   9.612,    0.150,  0.148,   0.152,    0.096,  0.086,   0.107,
  "Age",         "BlackBoost", 6.424,  3.333,   10.151,   0.574,  0.554,   0.587,    0.267,  0.247,   0.301,
  
  "Sex",         "Cox",        1.265,  0.346,   1.944,    2.194,  2.154,   2.245,    0.525,  0.474,   0.596,
  "Sex",         "RSF",        1.399,  0.410,   2.177,    0.150,  0.148,   0.152,    0.040,  0.036,   0.042,
  "Sex",         "BlackBoost", 1.454,  0.284,   2.460,    0.574,  0.554,   0.587,    0.133,  0.117,   0.151,
  
  "Income",      "Cox",        2.731,  0.563,   4.600,    2.194,  2.154,   2.245,    0.204,  0.138,   0.290,
  "Income",      "RSF",        2.999,  1.088,   4.633,    0.150,  0.148,   0.152,    0.024,  0.019,   0.029,
  "Income",      "BlackBoost", 2.611,  0.520,   4.619,    0.574,  0.554,   0.587,    0.086,  0.066,   0.108,
  
  "Marital",     "Cox",        0.721,  0.006,   2.416,    2.194,  2.154,   2.245,    0.399,  0.336,   0.438,
  "Marital",     "RSF",        0.783,  0.306,   1.452,    0.150,  0.148,   0.152,    0.030,  0.026,   0.036,
  "Marital",     "BlackBoost", 0.747,  0.057,   1.959,    0.574,  0.554,   0.587,    0.091,  0.086,   0.105,
  
  "Race",        "Cox",        1.599,  0.314,   3.000,    2.194,  2.154,   2.245,    1.002,  0.905,   1.054,
  "Race",        "RSF",        2.388,  1.084,   4.105,    0.150,  0.148,   0.152,    0.091,  0.086,   0.102,
  "Race",        "BlackBoost", 2.212,  0.279,   4.007,    0.574,  0.554,   0.587,    0.256,  0.228,   0.284,
  
  "Urban/Rural", "Cox",        1.137,  0.199,   2.272,    2.194,  2.154,   2.245,    0.097,  0.044,   0.150,
  "Urban/Rural", "RSF",        1.343,  0.213,   2.845,    0.150,  0.148,   0.152,    0.012,  0.006,   0.023,
  "Urban/Rural", "BlackBoost", 1.108,  0.141,   2.221,    0.574,  0.554,   0.587,    0.038,  0.008,   0.066
)

# Uzun formata getir (CI, Fi, Fg ayrı metric olacak)
ci_df <- fair_raw %>%
  select(Subgroup, Model, est = CI, lower = CI_low, upper = CI_high) %>%
  mutate(Metric = "CI%")

fi_df <- fair_raw %>%
  select(Subgroup, Model, est = Fi, lower = Fi_low, upper = Fi_high) %>%
  mutate(Metric = "Fi")

fg_df <- fair_raw %>%
  select(Subgroup, Model, est = Fg, lower = Fg_low, upper = Fg_high) %>%
  mutate(Metric = "Fg")

fair_long <- bind_rows(ci_df, fi_df, fg_df)

# Faktör sıraları (görüntü için)
fair_long <- fair_long %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c("Age", "Sex", "Income", "Marital", "Race", "Urban/Rural")),
    Model = factor(Model, levels = c("Cox", "RSF", "BlackBoost")),
    Metric = factor(Metric, levels = c("CI%", "Fi", "Fg"))
  )
library(ggplot2)

pd <- position_dodge(width = 0.7)

p_fair <- ggplot(fair_long,
                 aes(x = Subgroup, y = est, fill = Model)) +
  geom_col(position = pd, width = 0.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = pd, width = 0.25, linewidth = 0.3) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    scales = "free_y",
    labeller = labeller(
      Metric = c(
        "CI%" = "Concordance imparity (CI%)",
        "Fi"  = "Fi (individual fairness)",
        "Fg"  = "Fg (group fairness)"
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "Cox"        = "#4E79A7",  
      "RSF"        = "#F28E2B", 
      "BlackBoost" = "#59A14F"
    ),
    name = "Model"
  ) +
  labs(
    x = NULL,
    y = "Fairness metric value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.title.y       = element_text(size = 10)
  )

p_fair




# PNG (raster)
ggsave("Figure_fairness_barpanels.png",
       p_fair,
       width = 8, height = 3, dpi = 300)

# PDF (vector)
ggsave("Figure_fairness_barpanels.pdf",
       p_fair,
       width = 8, height = 3,
       device = "pdf")   



#### bar charts without errors ###
library(tidyverse)
library(ggplot2)

# ============================================================
# 1) Raw data (convert to tidy)
# ============================================================

fair_raw <- tribble(
  ~Subgroup,     ~Model,       ~CI,    ~CI_low, ~CI_high, ~Fi,    ~Fi_low, ~Fi_high, ~Fg,    ~Fg_low, ~Fg_high,
  "Age",         "CoxPH",        5.980,  3.649,   8.903,    2.194,  2.154,   2.245,    1.127,  0.956,   1.527,
  "Age",         "RSF",        5.946,  3.038,   9.612,    0.150,  0.148,   0.152,    0.096,  0.086,   0.107,
  "Age",         "BlackBoost", 6.424,  3.333,   10.151,   0.574,  0.554,   0.587,    0.267,  0.247,   0.301,
  
  "Sex",         "CoxPH",        1.265,  0.346,   1.944,    2.194,  2.154,   2.245,    0.525,  0.474,   0.596,
  "Sex",         "RSF",        1.399,  0.410,   2.177,    0.150,  0.148,   0.152,    0.040,  0.036,   0.042,
  "Sex",         "BlackBoost", 1.454,  0.284,   2.460,    0.574,  0.554,   0.587,    0.133,  0.117,   0.151,
  
  "Income",      "CoxPH",        2.731,  0.563,   4.600,    2.194,  2.154,   2.245,    0.204,  0.138,   0.290,
  "Income",      "RSF",        2.999,  1.088,   4.633,    0.150,  0.148,   0.152,    0.024,  0.019,   0.029,
  "Income",      "BlackBoost", 2.611,  0.520,   4.619,    0.574,  0.554,   0.587,    0.086,  0.066,   0.108,
  
  "Marital",     "CoxPH",        0.721,  0.006,   2.416,    2.194,  2.154,   2.245,    0.399,  0.336,   0.438,
  "Marital",     "RSF",        0.783,  0.306,   1.452,    0.150,  0.148,   0.152,    0.030,  0.026,   0.036,
  "Marital",     "BlackBoost", 0.747,  0.057,   1.959,    0.574,  0.554,   0.587,    0.091,  0.086,   0.105,
  
  "Race",        "CoxPH",        1.599,  0.314,   3.000,    2.194,  2.154,   2.245,    1.002,  0.905,   1.054,
  "Race",        "RSF",        2.388,  1.084,   4.105,    0.150,  0.148,   0.152,    0.091,  0.086,   0.102,
  "Race",        "BlackBoost", 2.212,  0.279,   4.007,    0.574,  0.554,   0.587,    0.256,  0.228,   0.284,
  
  "Urban/Rural", "CoxPH",        1.137,  0.199,   2.272,    2.194,  2.154,   2.245,    0.097,  0.044,   0.150,
  "Urban/Rural", "RSF",        1.343,  0.213,   2.845,    0.150,  0.148,   0.152,    0.012,  0.006,   0.023,
  "Urban/Rural", "BlackBoost", 1.108,  0.141,   2.221,    0.574,  0.554,   0.587,    0.038,  0.008,   0.066
)

# ---- reshape into long tidy form ----
ci_df <- fair_raw %>% 
  select(Subgroup, Model, est = CI, Metric = CI)

fi_df <- fair_raw %>% 
  select(Subgroup, Model, est = Fi) %>%
  mutate(Metric = "Fi")

fg_df <- fair_raw %>% 
  select(Subgroup, Model, est = Fg) %>%
  mutate(Metric = "Fg")

fair_long <- fair_raw %>%
  select(Subgroup, Model, CI, Fi, Fg) %>%
  pivot_longer(cols = c(CI, Fi, Fg), 
               names_to = "Metric", 
               values_to = "est")

# Factor ordering
fair_long <- fair_long %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c("Age", "Sex", "Income", "Marital", "Race", "Urban/Rural")),
    Model = factor(Model, 
                   levels = c("CoxPH", "RSF", "BlackBoost")),
    Metric = factor(Metric,
                    levels = c("CI", "Fi", "Fg"))
  )

# ============================================================
# 2) Nature-style COLORS
# ============================================================

nature_colors <- c(
  "CoxPH"        = "#4E79A7",  # muted blue
  "RSF"        = "#F28E2B",  # muted orange
  "BlackBoost" = "#59A14F"   # muted green
)

# ============================================================
# 3) Plot (NO error bars)
# ============================================================

pd <- position_dodge(width = 0.7)

p_fair <- ggplot(fair_long, aes(x = Subgroup, y = est, fill = Model)) +
  geom_col(position = pd, width = 0.6) +
  facet_wrap(
    ~ Metric,
    nrow = 1,
    scales = "free_y",
    labeller = labeller(
      Metric = c(
        "CI" = "Concordance imparity (CI%)",
        "Fi" = "Fi (individual fairness)",
        "Fg" = "Fg (group fairness)"
      )
    )
  ) +
  scale_fill_manual(values = nature_colors, name = "Model") +
  labs(x = NULL, y = "Fairness metric value") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "top",
    legend.title       = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

p_fair

# ============================================================
# 4) Save (PNG + Vector PDF)
# ============================================================

# PNG
ggsave("Figure_fairness_barpanels2.png",
       p_fair, width = 8, height = 3, dpi = 300)

# Vector PDF
ggsave("Figure_fairness_barpanels2.pdf",
       p_fair, width = 8, height = 3, device = "pdf")




##### Final code #####
##### Final Fairness Barplot Code #####
library(tidyverse)
library(ggplot2)

# ============================================================
# 1) Raw data
# ============================================================

fair_raw <- tribble(
  ~Subgroup,     ~Model,       ~CI,    ~CI_low, ~CI_high, ~Fi,    ~Fi_low, ~Fi_high, ~Fg,    ~Fg_low, ~Fg_high,
  
  "Age",         "CoxPH",        5.980,  3.649,   8.903,    2.194,  2.154,   2.245,    1.127,  0.956,   1.527,
  "Age",         "RSF",          5.946,  3.038,   9.612,    0.150,  0.148,   0.152,    0.096,  0.086,   0.107,
  "Age",         "BlackBoost",   6.424,  3.333,   10.151,   0.574,  0.554,   0.587,    0.267,  0.247,   0.301,
  
  "Sex",         "CoxPH",        1.265,  0.346,   1.944,    2.194,  2.154,   2.245,    0.525,  0.474,   0.596,
  "Sex",         "RSF",          1.399,  0.410,   2.177,    0.150,  0.148,   0.152,    0.040,  0.036,   0.042,
  "Sex",         "BlackBoost",   1.454,  0.284,   2.460,    0.574,  0.554,   0.587,    0.133,  0.117,   0.151,
  
  "Race",        "CoxPH",        1.599,  0.314,   3.000,    2.194,  2.154,   2.245,    1.002,  0.905,   1.054,
  "Race",        "RSF",          2.388,  1.084,   4.105,    0.150,  0.148,   0.152,    0.091,  0.086,   0.102,
  "Race",        "BlackBoost",   2.212,  0.279,   4.007,    0.574,  0.554,   0.587,    0.256,  0.228,   0.284,
  
  "Marital",     "CoxPH",        0.721,  0.006,   2.416,    2.194,  2.154,   2.245,    0.399,  0.336,   0.438,
  "Marital",     "RSF",          0.783,  0.306,   1.452,    0.150,  0.148,   0.152,    0.030,  0.026,   0.036,
  "Marital",     "BlackBoost",   0.747,  0.057,   1.959,    0.574,  0.554,   0.587,    0.091,  0.086,   0.105,
  
  "Income",      "CoxPH",        2.731,  0.563,   4.600,    2.194,  2.154,   2.245,    0.204,  0.138,   0.290,
  "Income",      "RSF",          2.999,  1.088,   4.633,    0.150,  0.148,   0.152,    0.024,  0.019,   0.029,
  "Income",      "BlackBoost",   2.611,  0.520,   4.619,    0.574,  0.554,   0.587,    0.086,  0.066,   0.108,
  
  "Urban/Rural", "CoxPH",        1.137,  0.199,   2.272,    2.194,  2.154,   2.245,    0.097,  0.044,   0.150,
  "Urban/Rural", "RSF",          1.343,  0.213,   2.845,    0.150,  0.148,   0.152,    0.012,  0.006,   0.023,
  "Urban/Rural", "BlackBoost",   1.108,  0.141,   2.221,    0.574,  0.554,   0.587,    0.038,  0.008,   0.066
)

# ============================================================
# 2) Build final long data
# ============================================================

# Subgroup metrics (CI, Fg)
ci_df <- fair_raw %>% transmute(Subgroup, Model, Metric = "CI", est = CI)
fg_df <- fair_raw %>% transmute(Subgroup, Model, Metric = "Fg", est = Fg)

# Fi: only one value per model
fi_df <- fair_raw %>%
  group_by(Model) %>%
  summarise(est = first(Fi), .groups = "drop") %>%
  mutate(Subgroup = "All patients",
         Metric   = "Fi") %>%
  select(Subgroup, Model, Metric, est)

# Combine all
fair_long <- bind_rows(fi_df, ci_df, fg_df)

# Correct ordering
fair_long <- fair_long %>%
  mutate(
    Subgroup = factor(Subgroup,
                      levels = c("All patients", "Age", "Sex", "Race", "Marital", "Income", "Urban/Rural")),
    Metric = factor(Metric, levels = c("Fi", "CI", "Fg")),
    Model = factor(Model, levels = c("CoxPH", "RSF", "BlackBoost"))
  )

# ============================================================
# 3) Plot
# ============================================================

nature_colors <- c("CoxPH"="#4E79A7","RSF"="#F28E2B","BlackBoost"="#59A14F")

p_fair <- ggplot(fair_long, aes(x = Subgroup, y = est, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ Metric, nrow = 1, scales = "free",
             labeller = labeller(
               Metric = c(
                 "Fi" = "* Fi (individual fairness)",
                 "CI" = "Concordance imparity (CI%)",
                 "Fg" = "Fg (group fairness)"
               )
             )
  ) +
  scale_fill_manual(values = nature_colors, name = "Model") +
  labs(x = NULL, y = "Fairness metric value") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_fair

# ============================================================
# 4) Save
# ============================================================

ggsave("Figure_fairness_bars_final3.png", p_fair,
       width = 8, height = 3, dpi = 300)

ggsave("Figure_fairness_bars_final3.pdf", p_fair,
       width = 8, height = 3, device = "pdf")
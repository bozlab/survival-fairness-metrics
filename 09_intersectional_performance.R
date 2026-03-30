# ============================================================
# FINAL SCRIPT — Intersectional cell-level C-index 
# Using saved RDS models/test only; no refitting.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, mc.cores = 1)  

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Hmisc)
  library(readr)
})

# -----------------------------
# 0) Paths (EDIT if needed)
# -----------------------------
models_rds <- "fitted_models_OS.rds"
test_rds   <- "test_data_seer_test.rds"
OUT_DIR    <- "figures_tauInf"

stopifnot(file.exists(models_rds), file.exists(test_rds))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Load saved objects
# -----------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

mods <- readRDS(models_rds)
test <- readRDS(test_rds)

cox_model <- mods$cox_model
rsf_model <- mods$rsf_model
bb_model  <- mods$bb_model
tau_test_inf <- mods$tau_test_inf %||% median(test$Survival.months, na.rm = TRUE)

# -----------------------------
# 2) Basic survival guards
# -----------------------------
complete_surv <- function(df){
  ok <- is.finite(df$Survival.months) & !is.na(df$event_status)
  df[ok, , drop = FALSE]
}

# -----------------------------
# 3) Risk scorer aligned with your pipeline
#    - Cox/BlackBoost: try risk/lp/link
#    - RSF: S3 predict.rfsrc; risk = 1 - S(t_med)
# -----------------------------
scalar_risk <- function(model, newdata){
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "risk")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "lp")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "link")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  
  if (inherits(model, "rfsrc")){
    med_t <- median(newdata$Survival.months, na.rm = TRUE)
    pr <- predict(model, newdata = newdata, na.action = "na.impute")  # predict.rfsrc
    j  <- which.min(abs(pr$time.interest - med_t))
    return(1 - as.numeric(pr$survival[, j]))
  }
  
  # Fallback (rare)
  r <- tryCatch({
    riskRegression::predictRisk(model, newdata = newdata, times = median(newdata$Survival.months, na.rm = TRUE))[,1]
  }, error = function(e) rep(NA_real_, nrow(newdata)))
  as.numeric(r)
}

# -----------------------------
# 4) Cell-level C-index at τ = ∞ (proxy via scalar_risk)
# -----------------------------
cell_cindex <- function(df, model, min_ok = 30){
  df <- complete_surv(df)
  if (nrow(df) < 50 || length(unique(df$event_status)) < 2) return(NA_real_)
  r <- scalar_risk(model, df)
  ok <- is.finite(r)
  if (sum(ok) < min_ok || length(unique(df$event_status[ok])) < 2) return(NA_real_)
  as.numeric(Hmisc::rcorr.cens(-r[ok], Surv(df$Survival.months[ok], df$event_status[ok]))["C Index"])
}

# -----------------------------
# 5) Settings / knobs
# -----------------------------
FAIR_VARS     <- c("Age_grouped","Race_grouped","Gender","Marital_grouped","Income_group","Urban_Rural")
MIN_CELL_N    <- 50              # lower to 25 if cells are sparse
RUN_ALL_PAIRS <- TRUE            # TRUE: all 2-way pairs; FALSE: only PAIRS_TO_RUN

PAIRS_TO_RUN  <- list(           # used only if RUN_ALL_PAIRS == FALSE
  c("Age_grouped","Income_group"),
  c("Race_grouped","Age_grouped"),
  c("Gender","Race_grouped")
)

# Display labels only (do not rename columns in data)
pretty_var <- function(v) ifelse(v == "Gender", "Sex", v)

MODELS_LIST   <- list(CoxPH = cox_model, RSF = rsf_model, BlackBoost = bb_model)

# -----------------------------
# 6) Compute one pair table (no bootstrap)
# -----------------------------
compute_pair <- function(var1, var2){
  if (!all(c(var1, var2) %in% names(test))) {
    message(sprintf("[skip] Missing columns: %s × %s", var1, var2)); return(NULL)
  }
  df <- test %>%
    filter(!is.na(.data[[var1]]), !is.na(.data[[var2]])) %>%
    mutate(
      !!var1 := factor(.data[[var1]]),
      !!var2 := factor(.data[[var2]])
    )
  if (nlevels(df[[var1]]) < 2 || nlevels(df[[var2]]) < 2) {
    message(sprintf("[skip] Not enough levels after filtering: %s × %s", var1, var2)); return(NULL)
  }
  
  grid <- expand.grid(L1 = levels(df[[var1]]), L2 = levels(df[[var2]]), stringsAsFactors = FALSE)
  
  rows <- lapply(seq_len(nrow(grid)), function(i){
    L1 <- grid$L1[i]; L2 <- grid$L2[i]
    sub <- df[df[[var1]] == L1 & df[[var2]] == L2, , drop = FALSE]
    n   <- nrow(sub)
    if (n < MIN_CELL_N) return(NULL)
    
    c_list <- lapply(MODELS_LIST, function(m) cell_cindex(sub, m))
    data.frame(
      Pair   = paste0(pretty_var(var1), " × ", pretty_var(var2)),
      Var1   = pretty_var(var1),
      Level1 = L1,
      Var2   = pretty_var(var2),
      Level2 = L2,
      n      = n,
      Model  = names(c_list),
      Cindex = as.numeric(c_list),
      row.names = NULL,
      check.names = FALSE
    )
  })
  
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) return(NULL)
  out
}

# -----------------------------
# 7) Run pairs
# -----------------------------
pair_list <- if (RUN_ALL_PAIRS) combn(FAIR_VARS, 2, simplify = FALSE) else PAIRS_TO_RUN
tables <- lapply(pair_list, function(p) compute_pair(p[1], p[2]))
inter_cells <- do.call(rbind, tables[!vapply(tables, is.null, logical(1))])

if (is.null(inter_cells) || !nrow(inter_cells)) {
  stop("No cells produced. Try lowering MIN_CELL_N or verify pair variable names.")
}

# -----------------------------
# 8) Worst cells per Pair & Model (lowest C-index)
# -----------------------------
worst_cells <- inter_cells %>%
  filter(is.finite(Cindex)) %>%
  group_by(Pair, Model) %>%
  slice_min(order_by = Cindex, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Pair, Model, Cindex)

cat("\n=== Worst 3 cells per Pair & Model (τ=∞) ===\n")
print(worst_cells %>% select(Pair, Model, Var1, Level1, Var2, Level2, n, Cindex),
      n = 200, row.names = FALSE)

# -----------------------------
# 9) Gap-to-best analysis
# -----------------------------
gap_tbl <- inter_cells %>%
  group_by(Pair, Model) %>%
  mutate(MaxC = suppressWarnings(max(Cindex, na.rm = TRUE)),
         MaxC = ifelse(is.finite(MaxC), MaxC, NA_real_),
         Gap  = ifelse(is.finite(Cindex) & is.finite(MaxC), MaxC - Cindex, NA_real_)) %>%
  ungroup()

largest_gaps <- gap_tbl %>%
  filter(is.finite(Gap)) %>%
  group_by(Pair, Model) %>%
  slice_max(order_by = Gap, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(Gap))

cat("\n=== Largest gap-to-best per Pair & Model (τ=∞) ===\n")
print(largest_gaps %>% select(Pair, Model, Var1, Level1, Var2, Level2, n, Cindex, MaxC, Gap),
      n = 200, row.names = FALSE)

# -----------------------------
# 10) Recurring vulnerable levels
# -----------------------------
vulnerable_levels <- worst_cells %>%
  mutate(Cell1 = paste0(Var1, "=", Level1),
         Cell2 = paste0(Var2, "=", Level2)) %>%
  pivot_longer(cols = c(Cell1, Cell2), names_to = "Axis", values_to = "LevelTag") %>%
  count(LevelTag, sort = TRUE)

cat("\n=== Levels most often appearing among worst cells (count across all pairs/models) ===\n")
print(vulnerable_levels, n = 50)

# -----------------------------
# 11) Heatmap utility + save all pairs
# -----------------------------
plot_cindex_heatmap <- function(tbl, pair_tag){
  pdat <- tbl %>% filter(Pair == pair_tag)
  if (!nrow(pdat)) return(NULL)
  pdat <- pdat %>% mutate(Cell = paste0(Var1, "=", Level1, "\n", Var2, "=", Level2))
  ggplot(pdat, aes(x = Model, y = Cell, fill = Cindex)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.finite(Cindex), sprintf("%.3f", Cindex), "NA")), size = 3) +
    scale_fill_gradient(limits = c(0.50, 0.90), low = "white", high = "steelblue", na.value = "grey90",
                        name = "C-index") +
    labs(title = paste0("Intersectional performance (τ=∞): ", pair_tag),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.title  = element_text(face = "bold"))
}

sanitize <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

unique_pairs <- unique(inter_cells$Pair)
for (pp in unique_pairs) {
  fig <- plot_cindex_heatmap(inter_cells, pp)
  if (!is.null(fig)) {
    out_png <- file.path(OUT_DIR, paste0("heatmap_cindex_", sanitize(pp), ".png"))
    ggsave(filename = out_png, plot = fig, width = 8, height = 7, dpi = 300)
    message("Saved figure: ", out_png)
  }
}

# -----------------------------
# 12) Save CSV outputs
# -----------------------------
write_csv(inter_cells,  file.path(OUT_DIR, "inter_cells_tauInf_ALL.csv"))
write_csv(worst_cells,  file.path(OUT_DIR, "inter_cells_tauInf_worst3.csv"))
write_csv(largest_gaps, file.path(OUT_DIR, "inter_cells_tauInf_largest_gaps.csv"))
write_csv(vulnerable_levels, file.path(OUT_DIR, "inter_cells_tauInf_vulnerable_levels.csv"))

cat("\n Export complete. Folder: ", OUT_DIR, "\n")
# ============================================================

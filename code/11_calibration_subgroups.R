# ============================================================
# Calibration-by-Subgroups (Style-matched) @ 24 / 36 / 60 months
#   * Loads models & test from RDS (no retraining)
#   * Normalizes labels, harmonizes factor levels for prediction
#   * ONE vector PDF per model (CoxPH / RSF / BlackBoost)
#       - 6 rows  = subgroup variables
#       - 3 cols  = time points (24, 36, 60 months)
#   * Axes locked to [0, AX_MAX], light grid, thin lines
#   * Robust binning + safe-time fallback so 60m never drops
# ============================================================

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(ggplot2)
  library(pec); library(randomForestSRC); library(patchwork); library(scales)
})

# --------------------------- Paths
MODELS_RDS <- "fitted_models_OS.rds"
TEST_RDS   <- "test_data_seer_test.rds"
PLOT_DIR   <- "plots_calibration"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------- What to plot
# Row order: Age → Sex → Race → Marital → Income → Urban/Rural
FAIR_VARS <- c("Age_grouped",
               "Gender",
               "Race_grouped",
               "Marital_grouped",
               "Income_group",
               "Urban_Rural")

ADD_OVERALL_ROW <- FALSE
AX_MAX          <- 0.80   # you can set to 1.0 if you want exactly the same range as 12-month plot
K_BINS          <- 10

# --------------------------- Desired legend order (pretty labels)
# NOTE: Race burada 12 ay kodundaki gibi bırakıldı (Black, Other, White)
#       Income ise 74,999 ve >$90,000 ile düzeltildi.
LEVEL_ORDER <- list(
  Age_grouped      = c("20–54 years","55–64 years","65–74 years","75–84 years","85+ years"),
  Income_group     = c("<$55,000","$55,000–$74,999","$75,000–$89,999",">$90,000"),
  Gender           = c("Female","Male"),
  Race_grouped     = c("Black","Other","White"),
  Marital_grouped  = c("Married","Others"),
  Urban_Rural      = c("Rural","Urban")
)

# Pretty display labels for legends (match 12-month figure)
DISPLAY_LABELS <- list(
  Gender          = "Sex",
  Race_grouped    = "Race",
  Age_grouped     = "Age",
  Marital_grouped = "Marital Status",
  Income_group    = "Income",
  Urban_Rural     = "Urban/Rural"
)

# If the model was trained with numeric-coded levels (e.g., "1","2",...),
# map codes <-> pretty labels here (both directions)
PRETTY_MAP <- list(
  Age_grouped_num      = c("1"="20–54 years","2"="55–64 years","3"="65–74 years","4"="75–84 years","5"="85+ years"),
  Income_group_num     = c("1"="<$55,000","2"="$55,000–$74,999","3"="$75,000–$89,999","4"=">$90,000"),
  Urban_Rural_num      = c("1"="Rural","2"="Urban"),
  Gender_num           = c("1"="Female","2"="Male"),
  Race_grouped_num     = c("1"="Black","2"="Other","3"="White"),
  Marital_grouped_num  = c("1"="Married","2"="Others")
)

# --------------------------- Palette & shapes
palette_levels <- function(n){
  base <- c("#E76F51","#2A9D8F","#F4A261","#457B9D","#9C27B0","#264653")
  base[seq_len(n)]
}
shape_levels <- function(n) rep(16, n)     # filled circles

# --------------------------- Load models & data
mods <- readRDS(MODELS_RDS)
cox_model <- mods$cox_model
rsf_model <- mods$rsf_model
bb_model  <- mods$bb_model
attr(bb_model, "bb_baseline") <- mods$bb_baseline

test_data <- readRDS(TEST_RDS)

# --------------------------- Normalize labels to canonical forms used in LEVEL_ORDER
normalize_labels <- function(df){
  trim_ws <- function(x) gsub("^\\s+|\\s+$", "", x)
  
  # Income_group: unify to canonical labels in LEVEL_ORDER
  if ("Income_group" %in% names(df)){
    inc <- as.character(df$Income_group)
    inc <- trim_ws(inc)
    inc <- gsub("<\\s*\\$55,000", "<$55,000", inc)
    inc <- gsub("\\$55,000\\s*-\\s*\\$75,000", "$55,000–$74,999", inc)
    inc <- gsub("\\$75,000\\s*-\\s*\\$89,999", "$75,000–$89,999", inc)
    inc <- gsub(">\\s*\\$90,000", ">$90,000", inc)
    df$Income_group <- factor(inc, levels = LEVEL_ORDER$Income_group)
  }
  
  # Age_grouped: ensure en-dash between ranges
  if ("Age_grouped" %in% names(df)){
    ag <- gsub(" ?- ?", "–", as.character(df$Age_grouped))
    df$Age_grouped <- factor(ag, levels = LEVEL_ORDER$Age_grouped)
  }
  
  # Fix explicit orders for the remaining factors
  if ("Urban_Rural" %in% names(df)){
    df$Urban_Rural <- factor(as.character(df$Urban_Rural), levels = LEVEL_ORDER$Urban_Rural)
  }
  if ("Gender" %in% names(df)){
    df$Gender <- factor(as.character(df$Gender), levels = LEVEL_ORDER$Gender)
  }
  if ("Race_grouped" %in% names(df)){
    df$Race_grouped <- factor(as.character(df$Race_grouped), levels = LEVEL_ORDER$Race_grouped)
  }
  if ("Marital_grouped" %in% names(df)){
    df$Marital_grouped <- factor(as.character(df$Marital_grouped), levels = LEVEL_ORDER$Marital_grouped)
  }
  df
}

# --------------------------- Survival completeness & tau
complete_surv <- function(df){
  ok <- is.finite(df$Survival.months) & !is.na(df$event_status)
  df[ok,,drop=FALSE]
}
choose_tau_inf <- function(df){
  tmax <- suppressWarnings(max(df$Survival.months, na.rm=TRUE))
  tevt <- suppressWarnings(max(df$Survival.months[df$event_status==1], na.rm=TRUE))
  tau <- min(tmax, tevt)
  if (!is.finite(tau) || tau <= 0) tau <- 36
  tau
}

# --------------------------- Harmonize to model xlevels for prediction
harmonize_to_model <- function(df, model){
  if (is.null(model$xlevels)) return(df)
  for (vn in intersect(names(model$xlevels), names(df))){
    x <- as.character(df[[vn]])
    x <- gsub(" ?- ?", "–", x, perl = TRUE)   # normalize dash to en-dash
    
    # Pretty -> numeric codes (when model$xlevels are "1","2",...)
    if (all(model$xlevels[[vn]] %in% c("1","2","3","4","5","6","7","8","9"))){
      if (vn == "Age_grouped"      && any(x %in% PRETTY_MAP$Age_grouped_num))     x <- names(PRETTY_MAP$Age_grouped_num)[match(x, PRETTY_MAP$Age_grouped_num)]
      if (vn == "Income_group"     && any(x %in% PRETTY_MAP$Income_group_num))    x <- names(PRETTY_MAP$Income_group_num)[match(x, PRETTY_MAP$Income_group_num)]
      if (vn == "Urban_Rural"      && any(x %in% PRETTY_MAP$Urban_Rural_num))     x <- names(PRETTY_MAP$Urban_Rural_num)[match(x, PRETTY_MAP$Urban_Rural_num)]
      if (vn == "Gender"           && any(x %in% PRETTY_MAP$Gender_num))          x <- names(PRETTY_MAP$Gender_num)[match(x, PRETTY_MAP$Gender_num)]
      if (vn == "Race_grouped"     && any(x %in% PRETTY_MAP$Race_grouped_num))    x <- names(PRETTY_MAP$Race_grouped_num)[match(x, PRETTY_MAP$Race_grouped_num)]
      if (vn == "Marital_grouped"  && any(x %in% PRETTY_MAP$Marital_grouped_num)) x <- names(PRETTY_MAP$Marital_grouped_num)[match(x, PRETTY_MAP$Marital_grouped_num)]
    }
    df[[vn]] <- factor(x, levels = model$xlevels[[vn]])
  }
  df
}

# --------------------------- RSF & BlackBoost adapters
predictSurvProb.blackboost <- function(object, newdata, times, ...){
  bl <- attr(object, "bb_baseline")
  if (is.null(bl) || is.null(bl$times) || is.null(bl$H0))
    stop("BlackBoost baseline not found (attach 'bb_baseline').")
  idx  <- findInterval(times, bl$times, left.open = FALSE)
  H0_t <- ifelse(idx == 0, 0, bl$H0[pmax(1, idx)])
  lp   <- as.numeric(predict(object, newdata = newdata, type = "link"))
  ee   <- exp(lp)
  S    <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(times))
  for (j in seq_along(times)) S[, j] <- exp(- H0_t[j] * ee)
  S
}
rsf_predict_surv_safe <- function(object, newdata, times){
  pr  <- predict(object, newdata = newdata, na.action = "na.impute")
  tt  <- pr$time.interest
  idx <- pmax(0L, findInterval(times, tt))
  S   <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  for (j in seq_along(times)) if (idx[j] > 0L) S[, j] <- pr$survival[, idx[j]]
  S
}

# --------------------------- Robust decile binning (handles flat predictions)
.make_bins <- function(S_vec, K=10){
  S_clean <- S_vec[is.finite(S_vec)]
  if (!length(S_clean)) return(rep(NA_integer_, length(S_vec)))
  br <- as.numeric(quantile(S_clean, probs = seq(0,1,length.out = K+1),
                            na.rm=TRUE, type=8))
  br <- unique(br)
  if (length(br) < 3){
    m <- min(S_clean, na.rm=TRUE); M <- max(S_clean, na.rm=TRUE)
    if (!is.finite(m) || !is.finite(M)) return(rep(NA_integer_, length(S_vec)))
    br <- c(m-1e-12, (m+M)/2, M+1e-12)
  } else {
    br <- cummax(br + 1e-12*seq_along(br))
  }
  cut(S_vec, breaks = br, include.lowest = TRUE, labels = FALSE)
}

# --------------------------- One subgroup level × one time t
calib_table_level <- function(df, S_vec, t, K=10){
  ok <- is.finite(S_vec) & is.finite(df$Survival.months) & !is.na(df$event_status)
  df <- df[ok,,drop=FALSE]; S_vec <- S_vec[ok]
  if (!nrow(df)) return(data.frame())
  bin <- .make_bins(S_vec, K); if (all(is.na(bin))) return(data.frame())
  
  km <- survfit(Surv(Survival.months, event_status) ~ bin, data=df)
  sm <- summary(km)  # full survival curve per bin
  
  nBins <- suppressWarnings(max(bin, na.rm=TRUE)); if (!is.finite(nBins)) return(data.frame())
  obs <- rep(NA_real_, nBins)
  for (b in seq_len(nBins)){
    idx <- which(sm$strata == paste0("bin=", b))
    if (!length(idx)) next
    tt <- sm$time[idx]; ss <- sm$surv[idx]
    t_eff <- min(t, max(tt, na.rm=TRUE))   # safe-time fallback
    if (all(tt > t_eff)) obs[b] <- 1.0 else obs[b] <- ss[max(which(tt <= t_eff))]
  }
  pred <- as.numeric(tapply(S_vec, bin, function(v) mean(v, na.rm=TRUE)))
  out  <- data.frame(Bin=seq_len(nBins), Predicted=pred, Observed=obs)
  out  <- out[is.finite(out$Predicted) & is.finite(out$Observed), , drop=FALSE]
  out[order(out$Predicted), , drop=FALSE]
}

# --------------------------- Build long table for all groups × times
build_long <- function(Smat, df, fair_vars, times, time_labels, K=10, add_overall=FALSE){
  pieces <- list()
  for (gvar in fair_vars){
    if (!gvar %in% names(df)) next
    levs <- levels(factor(df[[gvar]]))
    # respect desired legend order if present
    if (!is.null(LEVEL_ORDER[[gvar]])) levs <- intersect(LEVEL_ORDER[[gvar]], levs)
    for (j in seq_along(times)){
      t  <- times[j]
      tg <- time_labels[j]
      for (L in levs){
        idx <- which(df[[gvar]] == L); if (!length(idx)) next
        tab <- calib_table_level(df[idx,,drop=FALSE], Smat[idx, j], t, K)
        if (!nrow(tab)) next
        tab$GroupVar <- gvar
        tab$Level    <- L
        tab$TimeTag  <- tg
        pieces[[length(pieces)+1]] <- tab
      }
    }
  }
  if (!length(pieces)) return(data.frame())
  d <- bind_rows(pieces)
  d$GroupVar <- factor(d$GroupVar, levels = fair_vars)
  d
}

# --------------------------- After building, convert numeric codes -> pretty labels for plotting
prettify_levels_for_plot <- function(d){
  d$Level <- as.character(d$Level)
  
  map_char <- function(x, lut) {
    xc <- as.character(x)
    hit <- xc %in% names(lut)
    xc[hit] <- unname(lut[xc[hit]])
    xc
  }
  
  ok <- d$GroupVar == "Age_grouped"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Age_grouped_num)
  
  ok <- d$GroupVar == "Income_group"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Income_group_num)
  
  ok <- d$GroupVar == "Urban_Rural"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Urban_Rural_num)
  
  ok <- d$GroupVar == "Gender"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Gender_num)
  
  ok <- d$GroupVar == "Race_grouped"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Race_grouped_num)
  
  ok <- d$GroupVar == "Marital_grouped"
  if (any(ok)) d$Level[ok] <- map_char(d$Level[ok], PRETTY_MAP$Marital_grouped_num)
  
  d
}

# --------------------------- One-row plot (3 columns), legend on the right
row_plot <- function(drow, ax_max = 0.50){
  gname <- unique(as.character(drow$GroupVar))
  stopifnot(length(gname) == 1)
  
  # Display label for the legend title (e.g., "Sex" instead of "Gender")
  display_gname <- if (!is.null(DISPLAY_LABELS[[gname]])) DISPLAY_LABELS[[gname]] else gname
  
  # Factor to the desired order for this row
  desired <- if (!is.null(LEVEL_ORDER[[gname]])) {
    intersect(LEVEL_ORDER[[gname]], unique(drow$Level))
  } else {
    unique(drow$Level)
  }
  # 🔹 Race için sadece PLOT sırasını Black → White → Other yapıyoruz
  if (gname == "Race_grouped") {
    desired <- intersect(c("Black","White","Other"), desired)
  }
  drow$Level <- factor(drow$Level, levels = desired)
  
  lv  <- levels(drow$Level)
  pal <- setNames(palette_levels(length(lv)), lv)
  shp <- setNames(shape_levels(length(lv)),   lv)
  
  # Only 24 / 36 / 60 months
  drow$TimeTag <- factor(drow$TimeTag,
                         levels = c("24 Months","36 Months","60 Months"))
  drow$PanelLabel <- factor(
    paste0(display_gname, " - ", drow$TimeTag),
    levels = paste0(display_gname, " - ",
                    c("24 Months","36 Months","60 Months"))
  )
  
  ggplot(drow, aes(Predicted, Observed,
                   color = Level, shape = Level, group = Level)) +
    # Reference line y = x
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.5, color = "grey60") +
    
    # Calibration curve
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.0) +
    
    # Axes range
    coord_cartesian(xlim = c(0, ax_max), ylim = c(0, ax_max), expand = FALSE) +
    
    # X-axis: 0.1 increments
    scale_x_continuous(
      breaks = seq(0, ax_max, by = 0.1),
      labels = number_format(accuracy = 0.1)
    ) +
    # Y-axis: 0.1 increments
    scale_y_continuous(
      breaks = seq(0, ax_max, by = 0.1),
      labels = number_format(accuracy = 0.1)
    ) +
    
    # Manual color/shape settings
    scale_color_manual(values = pal, breaks = lv, limits = lv,
                       name = display_gname) +
    scale_shape_manual(values = shp,  breaks = lv, limits = lv,
                       name = display_gname) +
    
    # Facet across time points (3 columns)
    facet_wrap(~ PanelLabel, ncol = 3, drop = FALSE) +
    
    # Axis labels (match 12-month figure)
    labs(
      x = "Predicted Survival Probability",
      y = "Observed Survival Probability"
    ) +
    
    # Theme settings (spacing tuned to avoid overlapping tick labels)
    theme_bw() +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      
      strip.text = element_text(size = 9, face = "bold"),
      
      axis.title = element_text(size = 9),
      axis.text  = element_text(
        size = 8,
        margin = margin(t = 2, r = 2, b = 2, l = 2)
      ),
      
      # Wider spacing between columns & rows to prevent axis overlap
      panel.spacing.x = unit(2,  "lines"),
      panel.spacing.y = unit(1.2,"lines"),
      
      legend.position = "right",
      legend.title    = element_text(size = 9, face = "bold"),
      legend.text     = element_text(size = 8),
      
      # No extra main title / caption (you will describe it in the figure legend)
      plot.margin     = margin(4, 12, 4, 4)
    )
}

# --------------------------- Export: one PDF page for each model
export_onepage <- function(model_name, dlong, outfile, ax_max = 0.50){
  if (!nrow(dlong)) {
    message("No data to plot for ", model_name)
    return(invisible(NULL))
  }
  dlong$GroupVar <- factor(dlong$GroupVar, levels = FAIR_VARS)
  
  rows <- lapply(levels(dlong$GroupVar), function(g){
    drow <- dlong %>% filter(GroupVar == g)
    if (nrow(drow) == 0) return(NULL)
    row_plot(drow, ax_max = ax_max)
  })
  rows <- Filter(Negate(is.null), rows)
  
  pdf(outfile, width = 11, height = 16,
      useDingbats = FALSE, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  
  # Just stack the 6 rows; no global title (figure caption will be in the paper)
  print(
    wrap_plots(rows, ncol = 1)
  )
  
  message("Saved: ", outfile)
}

# --------------------------- Prep & predict

# 1) Keep only rows with time + event available
test_data <- complete_surv(test_data)

# 2) Normalize labels so Income/Gender etc. match LEVEL_ORDER canonically
test_data <- normalize_labels(test_data)

# 3) Harmonize to each model BEFORE predicting (handles numeric-coded levels safely)
df_cox <- harmonize_to_model(test_data, cox_model)
df_rsf <- harmonize_to_model(test_data, rsf_model)
df_bb  <- harmonize_to_model(test_data, bb_model)

# 4) Times (ONLY 24 / 36 / 60 months; no tau-infinity panel)
tau_inf   <- choose_tau_inf(test_data)
times_use <- c(min(24, tau_inf), min(36, tau_inf), min(60, tau_inf))
time_tags <- c("24 Months","36 Months","60 Months")
message(sprintf("τ∞ = %.2f | times = %s",
                tau_inf, paste(times_use, collapse = ", ")))

# 5) Predict survival probabilities
S_cox <- as.matrix(pec::predictSurvProb(cox_model, newdata = df_cox, times = times_use))
S_rsf <- as.matrix(rsf_predict_surv_safe(rsf_model, newdata = df_rsf, times = times_use))
S_bb  <- as.matrix(predictSurvProb(bb_model,   newdata = df_bb,  times = times_use))

# 6) Build long frames (+ prettify labels)
make_long <- function(S, df)
  build_long(S, df, FAIR_VARS, times_use, time_tags,
             K = K_BINS, add_overall = ADD_OVERALL_ROW)

d_cox <- prettify_levels_for_plot(make_long(S_cox, df_cox))
d_rsf <- prettify_levels_for_plot(make_long(S_rsf, df_rsf))
d_bb  <- prettify_levels_for_plot(make_long(S_bb,  df_bb))

# 7) Export PDFs (vector, one per model)
export_onepage("CoxPH",
               d_cox,
               file.path(PLOT_DIR, "Calibration_BySubgroups_CoxPH_styleMatched.pdf"),
               ax_max = AX_MAX)

export_onepage("RSF",
               d_rsf,
               file.path(PLOT_DIR, "Calibration_BySubgroups_RSF_styleMatched.pdf"),
               ax_max = AX_MAX)

export_onepage("BlackBoost",
               d_bb,
               file.path(PLOT_DIR, "Calibration_BySubgroups_BlackBoost_styleMatched.pdf"),
               ax_max = AX_MAX)

message("Saved 24/36/60-month calibration PDFs to: ", PLOT_DIR)
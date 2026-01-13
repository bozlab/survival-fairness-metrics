# =============================================================================
# Survival Fairness (Keya et al., 2021) — FULL FINAL SCRIPT (Train + Fairness)
#   - Trains CoxPH, RSF, BlackBoost from CSV (same seeds for reproducibility)
#   - Computes Keya 3.1 fairness metrics with bootstrap 95% CIs:
#       * Individual fairness  (Fi)     -> computed ONCE per sample (not per group)
#       * Group fairness       (Fg)
#       * Intersectional F∩    (2-way intersections)
#   - Outputs a single Excel workbook with safe sheet names.
#
#   QUICK RUN:   BOOT_N = 200,  RSF ntree = 500
#   FINAL RUN:   BOOT_N = 1000, RSF ntree = 1000  (keep seeds!)
#
#   Notes:
#   * Fi is group-agnostic by definition. We compute it once per (boot) sample
#     and replicate the same value across Sensitive to keep the wide-table layout.
#   * Determinism/stability: L'Ecuyer-CMRG RNG, mc.reset.stream(), mc.set.seed=TRUE,
#     and a larger Fi pair sample cap (pair_cap=100000L).
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(Hmisc)
  library(mboost)
  library(randomForestSRC)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(writexl)
  library(parallel)
  library(broom)
})

# --------------------------- Paths (EDIT to your actual CSV locations)
BASE       <- normalizePath("~/Library/CloudStorage/Dropbox/MyProject/TLCR", mustWork = TRUE)
TRAIN_CSV  <- file.path(BASE, "train_data_seer.csv")
TEST_CSV   <- file.path(BASE, "test_data_seer.csv")
OUT_XLS    <- file.path(BASE, "Fairness_Keya_Final3.xlsx")

# --------------------------- Settings (stable RNG & conservative parallelism)
set.seed(123); RNGkind("L'Ecuyer-CMRG")   # deterministic parallel RNG
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
if (.Platform$OS.type == "unix") parallel::mc.reset.stream()  # stable streams across forks

NCORES   <- min(6L, max(1L, parallel::detectCores() - 2L))
ALPHA    <- 0.05
TOP_K    <- 5
BOOT_N   <- 1000        # final paper run (use 50/200 for smoke tests)
C_CONST  <- 0.01        # Keya scaling constant (gamma)
MIN_N    <- 25          # minimum subgroup size for intersectional cells
FAIR_VARS <- c("Gender","Race_grouped","Age_grouped","Marital_grouped","Income_group","Urban_Rural")

# RSF settings (final: ntree = 1000)
RSF_NTREE    <- 1000
RSF_NODESIZE <- 15

# --------------------------- Small helpers

# 3-decimal string formatter (prints "NA" if non-finite)
fmt3 <- function(x) ifelse(is.finite(x), sprintf("%.3f", round(as.numeric(x), 3)), "NA")

# Convert to {0,1} with NA -> 0
to01 <- function(x){
  if (is.factor(x)) x <- as.character(x)
  x <- as.integer(x); x[is.na(x)] <- 0L; x[x != 0L] <- 1L; x
}

# Keep rows with finite survival time and event_status
complete_surv <- function(df){
  ok <- is.finite(df$Survival.months) & is.finite(df$event_status)
  df[ok,,drop=FALSE]
}

# --------------------------- Survival predictions (robust)

# RSF survival prediction on a requested time grid, preserving row order
rsf_predict_surv_safe <- function(object, newdata, times){
  pr <- predict(object, newdata = newdata, na.action = "na.impute")
  tt <- pr$time.interest
  idx <- pmax(0L, findInterval(times, tt))
  S <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  for (j in seq_along(times)) if (idx[j] > 0L) S[, j] <- pr$survival[, idx[j]]
  S
}

# Compute and attach BlackBoost baseline (Nelson–Aalen) for survival curves
.compute_bb_baseline <- function(model, time, status, newdata_for_eta){
  eta <- as.numeric(predict(model, newdata = newdata_for_eta, type = "link"))
  ord <- order(time, na.last = NA); time <- time[ord]; status <- status[ord]; eta <- eta[ord]
  ev_times <- sort(unique(time[status == 1])); H0 <- numeric(length(ev_times)); cum <- 0
  for (k in seq_along(ev_times)){
    tk <- ev_times[k]; d_k <- sum(status == 1 & time == tk); rset <- sum(exp(eta[time >= tk]))
    if (rset > 0) cum <- cum + d_k / rset; H0[k] <- cum
  }
  list(times = ev_times, H0 = H0)
}

predictSurvProb.blackboost <- function(object, newdata, times, ...){
  bl <- attr(object, "bb_baseline")
  if (is.null(bl) || is.null(bl$times) || is.null(bl$H0)) stop("BlackBoost baseline not found.")
  idx  <- findInterval(times, bl$times, left.open = FALSE)
  H0_t <- ifelse(idx == 0, 0, bl$H0[pmax(1, idx)])
  lp   <- as.numeric(predict(object, newdata = newdata, type = "link"))
  exp_eta <- exp(lp)
  N <- nrow(newdata); TT <- length(times); S <- matrix(NA_real_, nrow = N, ncol = TT)
  for (j in seq_len(TT)) S[, j] <- exp(- H0_t[j] * exp_eta)
  S
}
predictSurvProb.mboost <- function(object, newdata, times, ...) predictSurvProb.blackboost(object, newdata, times, ...)

# --------------------------- Unified risk proxy h(x) per Keya 3.1
# Cox/BlackBoost: h(x) = exp(linear predictor)
# RSF:            h(x) = mean_t [1 - S(t|x)] over a time grid
hx_score <- function(model, df, times_vec) {
  if (inherits(model, "coxph")) {
    lp <- as.numeric(predict(model, newdata = df, type = "lp"))
    return(exp(lp))
  }
  if (inherits(model, "mboost")) {
    lp <- as.numeric(predict(model, newdata = df, type = "link"))
    return(exp(lp))
  }
  if (inherits(model, "rfsrc")) {
    S <- rsf_predict_surv_safe(model, df, times_vec)
    return(drop(rowMeans(1 - S, na.rm = TRUE)))
  }
  stop("Unsupported model class in hx_score()")
}

# --------------------------- Row-preserving, NA-safe encoder (no contrasts)
# Creates an Euclidean design matrix on encoded features and scales it.
design_matrix_std <- function(df, drop_cols = c("Survival.months","event_status")) {
  keep <- setdiff(names(df), drop_cols)
  n <- nrow(df)
  if (!length(keep) || n == 0) return(matrix(0, nrow = n, ncol = 1))
  
  df2 <- df[, keep, drop = FALSE]
  
  # Safe column-wise scaler
  .scale_safe <- function(M){
    if (!ncol(M)) return(matrix(0, nrow = nrow(M), ncol = 1))
    mn  <- suppressWarnings(colMeans(M))
    sdv <- suppressWarnings(apply(M, 2, sd))
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    M <- sweep(M, 2, mn, FUN = "-")
    M <- sweep(M, 2, sdv, FUN = "/")
    M
  }
  
  pieces <- list()
  
  for (nm in colnames(df2)) {
    x <- df2[[nm]]
    if (all(is.na(x))) next  # skip all-NA columns
    
    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      v <- as.numeric(x)
      v[!is.finite(v)] <- 0  # numeric NA -> 0 before scaling
      pieces[[length(pieces)+1]] <- matrix(v, nrow = n, ncol = 1,
                                           dimnames = list(NULL, nm))
    } else {
      # character/factor: make NA an explicit level, then manual one-hot (no contrasts)
      f <- if (is.factor(x)) x else factor(x)
      f <- addNA(f, ifany = TRUE)
      levs <- levels(f); k <- length(levs)
      Z <- matrix(0, nrow = n, ncol = k)
      idx <- as.integer(f); ok <- is.finite(idx)
      Z[cbind(which(ok), idx[ok])] <- 1
      colnames(Z) <- paste0(nm, "=", levs)
      pieces[[length(pieces)+1]] <- Z
    }
  }
  
  if (!length(pieces)) return(matrix(0, nrow = n, ncol = 1))
  mm <- do.call(cbind, pieces)
  
  # Drop zero-variance columns
  if (ncol(mm) > 1) {
    vars <- apply(mm, 2, function(col) var(col, na.rm = TRUE))
    keep_cols <- is.finite(vars) & vars > 0
    if (any(!keep_cols)) mm <- mm[, keep_cols, drop = FALSE]
  }
  if (!ncol(mm)) mm <- matrix(0, nrow = n, ncol = 1)
  
  mm <- .scale_safe(mm)
  
  # Final safeguard: enforce row count equality
  if (nrow(mm) != n) {
    M2 <- matrix(0, nrow = n, ncol = ncol(mm))
    r  <- min(nrow(mm), n)
    M2[seq_len(r), ] <- mm[seq_len(r), , drop = FALSE]
    mm <- M2
  }
  mm
}

# --------------------------- Keya 3.1 metrics (point estimators)

# Individual fairness: Fi = mean_{i<j} max(0, |h_i - h_j| - C * ||x_i - x_j||_2)
# Sampling is capped for speed; index sampling is vectorized and i<j guaranteed.
Fi_point <- function(h, Xstd, C = C_CONST, pair_cap = 100000L) {  # ↑ larger cap for stability
  n <- length(h); if (n < 2) return(NA_real_)
  stopifnot(nrow(Xstd) == n)
  total_pairs <- n * (n - 1L) / 2L
  
  if (total_pairs > pair_cap) {
    # i ~ Uniform{1..n-1}
    i_idx <- sample.int(n - 1L, pair_cap, replace = TRUE)
    # j ~ Uniform{i+1 .. n}  --> offset in {1 .. n - i}
    offs  <- 1L + floor(runif(pair_cap) * (n - i_idx))
    j_idx <- i_idx + offs
  } else {
    ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
    i_idx <- ij[, 1]; j_idx <- ij[, 2]
  }
  if (!all(j_idx > i_idx & j_idx <= n & i_idx >= 1L)) stop("Fi_point: invalid index sampling.")
  
  dh <- abs(h[i_idx] - h[j_idx])
  dX <- sqrt(rowSums((Xstd[i_idx, , drop = FALSE] - Xstd[j_idx, , drop = FALSE])^2))
  mean(pmax(dh - C * dX, 0), na.rm = TRUE)
}

# Group fairness: Fg = max_a | E[h | A=a] - E[h] |
Fg_point <- function(h, A) {
  A <- factor(A)
  grp_means <- tapply(h, A, function(v) mean(v, na.rm = TRUE))
  pop_mean  <- mean(h, na.rm = TRUE)
  max(abs(grp_means - pop_mean), na.rm = TRUE)
}

# Intersectional fairness: F∩ = max_{si,sj} | log E[h | si] - log E[h | sj] |
Fcap_point <- function(h, df, vars, min_n = MIN_N) {
  keep <- stats::complete.cases(df[, vars, drop = FALSE])
  if (!any(keep)) return(NA_real_)
  df2 <- df[keep,,drop=FALSE]; h2 <- h[keep]
  key <- interaction(df2[, vars, drop = FALSE], drop = TRUE, sep = "×")
  tab <- split(h2, key)
  means <- vapply(tab, function(v) if (length(v) >= min_n) mean(v, na.rm = TRUE) else NA_real_, numeric(1))
  means <- means[is.finite(means) & means > 0]
  if (length(means) < 2) return(NA_real_)
  max(abs(outer(log(means), log(means), `-`)), na.rm = TRUE)
}

build_intersections <- function(vars) combn(vars, 2, simplify = FALSE)

# --------------------------- Load data
stopifnot(file.exists(TRAIN_CSV), file.exists(TEST_CSV))
train_data <- read.csv(TRAIN_CSV, check.names = TRUE)
test_data  <- read.csv(TEST_CSV,  check.names = TRUE)

# Outcomes
if (!"event_status" %in% names(train_data) && "Vital.status.recode" %in% names(train_data))
  train_data$event_status <- ifelse(train_data$Vital.status.recode=="Dead",1,0)
if (!"event_status" %in% names(test_data) && "Vital.status.recode" %in% names(test_data))
  test_data$event_status <- ifelse(test_data$Vital.status.recode=="Dead",1,0)

train_data$event_status    <- to01(train_data$event_status)
test_data$event_status     <- to01(test_data$event_status)
train_data$Survival.months <- as.numeric(train_data$Survival.months)
test_data$Survival.months  <- as.numeric(test_data$Survival.months)

# --------------------------- Predictors (exclude Income_group & Urban_Rural from model)
predictors <- c(
  "Age_grouped","Race_grouped","Gender","Marital_grouped",
  "Primary_site_grouped","Histology_grouped","Grade_grouped","Laterality_grouped",
  "T_stage_recoded","N_stage_combined","M_stage_combined",
  "Surgery","Chemotherapy","Radiation",
  "Bone_met","Liver_met","Lung_met","Brain_met"
)

# Coerce to factors
for (v in predictors){
  if (!is.null(train_data[[v]]) && !is.numeric(train_data[[v]])) train_data[[v]] <- factor(train_data[[v]])
  if (!is.null(test_data[[v]])  && !is.numeric(test_data[[v]]))  test_data[[v]]  <- factor(test_data[[v]])
}
train_data <- droplevels(train_data); test_data <- droplevels(test_data)

# Align factor levels across train/test for selected predictors
align_factors <- function(train, test, vars){
  for (v in vars){
    if (is.character(train[[v]])) train[[v]] <- factor(train[[v]])
    if (is.character(test[[v]]))  test[[v]]  <- factor(test[[v]])
    if (is.factor(train[[v]]) || is.factor(test[[v]])){
      lv <- union(levels(factor(train[[v]])), levels(factor(test[[v]])))
      train[[v]] <- factor(train[[v]], levels = lv)
      test[[v]]  <- factor(test[[v]],  levels = lv)
    }
  }
  list(train = train, test = test)
}

# --------------------------- Univariate screen to pick signals (same design)
surv_obj_tr <- with(train_data, Surv(Survival.months, event_status))
uni_results <- lapply(predictors, function(var){
  fit <- coxph(as.formula(paste("surv_obj_tr ~", sprintf("`%s`", var))), data = train_data)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(variable = var)
})
uni_df <- bind_rows(uni_results)
var_min_p <- uni_df %>% group_by(variable) %>% summarise(min_p = suppressWarnings(min(p.value, na.rm = TRUE)), .groups = "drop")

sig_vars <- var_min_p %>% filter(is.finite(min_p) & min_p < ALPHA) %>% pull(variable)
if (!length(sig_vars)){
  sig_vars <- var_min_p %>% arrange(min_p) %>% slice_head(n = min(TOP_K, n())) %>% pull(variable)
}
al <- align_factors(train_data, test_data, sig_vars); train_data <- al$train; test_data <- al$test

cox_formula <- as.formula(paste("Surv(Survival.months, event_status) ~", paste(sprintf("`%s`", sig_vars), collapse = " + ")))
message("[INFO] Cox formula: ", deparse(cox_formula))

# --------------------------- Fit models
cox_model <- coxph(cox_formula, data = train_data, x = TRUE, y = TRUE, model = TRUE)

p <- length(sig_vars)
rsf_model <- randomForestSRC::rfsrc(
  cox_formula, data = train_data,
  ntree = RSF_NTREE, splitrule = "logrank",
  mtry = max(1, floor(sqrt(p))), nodesize = RSF_NODESIZE, nsplit = 0,
  importance = "none", na.action = "na.impute"
)

bb_model <- mboost::blackboost(cox_formula, data = train_data, family = mboost::CoxPH())
set.seed(123)
cvb <- mboost::cvrisk(bb_model, folds = mboost::cv(model.weights(bb_model), type = "kfold", B = 5))
mboost::mstop(bb_model) <- mboost::mstop(cvb)
bb_baseline <- .compute_bb_baseline(bb_model, train_data$Survival.months, train_data$event_status, train_data)
attr(bb_model, "bb_baseline") <- bb_baseline

# --------------------------- Time grid for RSF/prob scores (from TEST set)
dfc <- complete_surv(test_data)
TIME_PROBS <- c(0.25, 0.50, 0.75)
tq <- as.numeric(quantile(dfc$Survival.months, probs = TIME_PROBS, na.rm = TRUE))
tq <- sort(unique(tq[is.finite(tq) & tq > 0])); if (!length(tq)) tq <- c(12, 24, 36)

# --------------------------- Bootstrap engine for Keya metrics (Fi once/replicate)
boot_keya <- function(df, model, model_name, fair_vars, times_vec, B = BOOT_N, ncores = NCORES) {
  df <- complete_surv(df)
  if (nrow(df) < 50) stop("Too few rows in test data after cleaning.")
  Xstd <- design_matrix_std(df)
  h    <- hx_score(model, df, times_vec)
  if (length(h) != nrow(Xstd)) {
    warning(sprintf("h (%d) vs Xstd rows (%d) mismatch — re-encoding.", length(h), nrow(Xstd)))
    Xstd <- design_matrix_std(df)
    if (length(h) != nrow(Xstd)) stop(sprintf("Still mismatched: h=%d, Xstd rows=%d", length(h), nrow(Xstd)))
  }
  # --- Fi ONCE on original sample ---
  Fi_pt_val <- Fi_point(h, Xstd, C = C_CONST)
  # --- Fg points (per sensitive); not shown in Excel but kept for traceability ---
  Fg_pt <- list()
  for (A in fair_vars) if (A %in% names(df)) Fg_pt[[A]] <- Fg_point(h, df[[A]])
  # --- Intersections to evaluate ---
  combos <- build_intersections(intersect(fair_vars, names(df)))
  Fcap_pt <- list()
  for (cmb in combos) { tag <- paste(cmb, collapse = "×"); Fcap_pt[[tag]] <- Fcap_point(h, df, cmb, min_n = MIN_N) }
  
  # --- Bootstrap worker (Fi ONCE per replicate) ---
  boot_core <- function(b) {
    idx  <- sample.int(nrow(df), replace = TRUE)
    boot <- complete_surv(df[idx,,drop=FALSE])
    if (nrow(boot) < 50) return(NULL)
    Xb <- design_matrix_std(boot)
    hb <- hx_score(model, boot, times_vec)
    if (length(hb) != nrow(Xb)) return(NULL)
    Fi_b_val <- Fi_point(hb, Xb, C = C_CONST)  # single scalar per replicate
    
    Fg_b <- list()
    for (A in fair_vars) if (A %in% names(boot)) Fg_b[[A]] <- Fg_point(hb, boot[[A]])
    Fcap_b <- list()
    for (cmb in combos) { tag <- paste(cmb, collapse = "×"); Fcap_b[[tag]] <- Fcap_point(hb, boot, cmb, min_n = MIN_N) }
    list(Fi = Fi_b_val, Fg = Fg_b, Fcap = Fcap_b)
  }
  
  # --- Chunked parallelization (deterministic seeding) ---
  chunks  <- split(seq_len(B), ceiling(seq_len(B) / max(1L, floor(B/100))))
  res_all <- list()
  for (ch in chunks) {
    res_ch <- parallel::mclapply(ch, boot_core, mc.cores = ncores, mc.set.seed = TRUE)
    res_all <- c(res_all, res_ch)
  }
  
  # --- Aggregate bootstrap stats (mean + percentile CI) ---
  agg_stat <- function(v) {
    v <- unlist(v, use.names = FALSE); v <- v[is.finite(v)]
    if (!length(v)) return(c(mean=NA, lower=NA, upper=NA, n=0))
    c(mean = mean(v), lower = as.numeric(quantile(v, 0.025)),
      upper = as.numeric(quantile(v, 0.975)), n = length(v))
  }
  
  # Fi as scalar -> replicate across Sensitive keys to keep layout
  sens_keys <- intersect(fair_vars, names(df)); if (!length(sens_keys)) sens_keys <- NA_character_
  fi_vals <- as.numeric(unlist(lapply(res_all, function(z) if(!is.null(z)) z$Fi)))
  fi_vals <- fi_vals[is.finite(fi_vals)]
  fi_stat <- if (!length(fi_vals)) c(mean=NA, lower=NA, upper=NA) else
    c(mean = mean(fi_vals), lower = as.numeric(quantile(fi_vals, 0.025)), upper = as.numeric(quantile(fi_vals, 0.975)))
  fi_tab <- data.frame(
    Method    = model_name,
    Sensitive = sens_keys,
    Fi_point  = Fi_pt_val,
    Fi_mean   = fi_stat["mean"],
    Fi_lower  = fi_stat["lower"],
    Fi_upper  = fi_stat["upper"],
    check.names = FALSE
  )
  
  # Fg per sensitive
  fg_keys <- sens_keys
  fg_tab <- bind_rows(lapply(fg_keys, function(A){
    st <- agg_stat(lapply(res_all, function(z) if(!is.null(z)) z$Fg[[A]] else NULL))
    data.frame(Method = model_name, Sensitive = A,
               Fg_point = Fg_pt[[A]], Fg_mean = st["mean"], Fg_lower = st["lower"], Fg_upper = st["upper"])
  }))
  
  # Fcap per intersection combo
  fc_keys <- vapply(combos, function(x) paste(x, collapse = "×"), character(1))
  fcap_tab <- bind_rows(lapply(fc_keys, function(K){
    st <- agg_stat(lapply(res_all, function(z) if(!is.null(z)) z$Fcap[[K]] else NULL))
    data.frame(Method = model_name, Combo = K,
               Fcap_point = Fcap_pt[[K]], Fcap_mean = st["mean"], Fcap_lower = st["lower"], Fcap_upper = st["upper"])
  }))
  
  list(Fi = fi_tab, Fg = fg_tab, Fcap = fcap_tab)
}

# --------------------------- Run for the three models
all_Fi <- list(); all_Fg <- list(); all_Fcap <- list()

message("[RUN] Keya fairness — Cox")
out <- boot_keya(test_data, cox_model, "Cox", FAIR_VARS, tq, B = BOOT_N, ncores = NCORES)
all_Fi[["Cox"]] <- out$Fi; all_Fg[["Cox"]] <- out$Fg; all_Fcap[["Cox"]] <- out$Fcap

message("[RUN] Keya fairness — RSF")
out <- boot_keya(test_data, rsf_model, "RSF", FAIR_VARS, tq, B = BOOT_N, ncores = NCORES)
all_Fi[["RSF"]] <- out$Fi; all_Fg[["RSF"]] <- out$Fg; all_Fcap[["RSF"]] <- out$Fcap

message("[RUN] Keya fairness — BlackBoost")
out <- boot_keya(test_data, bb_model, "BlackBoost", FAIR_VARS, tq, B = BOOT_N, ncores = NCORES)
all_Fi[["BlackBoost"]] <- out$Fi; all_Fg[["BlackBoost"]] <- out$Fg; all_Fcap[["BlackBoost"]] <- out$Fcap

Fi_tab   <- bind_rows(all_Fi)
Fg_tab   <- bind_rows(all_Fg)
Fcap_tab <- bind_rows(all_Fcap)

# --------------------------- Format (bootstrap-centric; no separate "point" column)
Fi_out <- Fi_tab %>%
  transmute(Method, Sensitive,
            `Fi [95% CI]` = paste0(fmt3(Fi_mean), " [", fmt3(Fi_lower), ", ", fmt3(Fi_upper), "]"))

Fg_out <- Fg_tab %>%
  transmute(Method, Sensitive,
            `Fg [95% CI]` = paste0(fmt3(Fg_mean), " [", fmt3(Fg_lower), ", ", fmt3(Fg_upper), "]"))

Fcap_out <- Fcap_tab %>%
  transmute(Method, Combo,
            `F∩ [95% CI]` = paste0(fmt3(Fcap_mean), " [", fmt3(Fcap_lower), ", ", fmt3(Fcap_upper), "]"))

# --------------------------- Detailed E[h | subgroup] (descriptive means; no CI)
intersection_means <- function(df, model, model_name, fair_vars, times_vec, min_n = MIN_N){
  df <- complete_surv(df); h <- hx_score(model, df, times_vec)
  combos <- build_intersections(intersect(fair_vars, names(df)))
  rows <- list()
  for (cmb in combos) {
    tag <- paste(cmb, collapse = "×")
    key <- interaction(df[, cmb, drop = FALSE], drop = TRUE, sep = "×")
    tab <- split(h, key)
    for (nm in names(tab)) {
      v <- tab[[nm]]
      if (length(v) >= min_n) {
        rows[[length(rows)+1]] <- data.frame(
          Method = model_name, Combo = tag, Subgroup = nm, n = length(v),
          `E[h|subgroup]` = mean(v, na.rm = TRUE), check.names = FALSE
        )
      }
    }
  }
  if (length(rows)) bind_rows(rows) else
    data.frame(Method=character(0), Combo=character(0), Subgroup=character(0), n=integer(0), `E[h|subgroup]`=numeric(0))
}

detail_list <- list(
  Cox        = intersection_means(test_data, cox_model,        "Cox",        FAIR_VARS, tq, MIN_N),
  RSF        = intersection_means(test_data, rsf_model,        "RSF",        FAIR_VARS, tq, MIN_N),
  BlackBoost = intersection_means(test_data, bb_model, "BlackBoost", FAIR_VARS, tq, MIN_N)
)
Detail_subgroups <- bind_rows(detail_list) %>%
  mutate(`E[h|subgroup] (3dp)` = fmt3(`E[h|subgroup]`)) %>%
  select(Method, Combo, Subgroup, n, `E[h|subgroup] (3dp)`)

# --------------------------- Write Excel with safe sheet names
sheets <- list(
  Fi_boot_CI   = Fi_out,
  Fg_boot_CI   = Fg_out,
  Fcap_boot_CI = Fcap_out,         # 2-way intersections; bootstrap mean + 95% CI
  Eh_by_subgrp = Detail_subgroups  # descriptive subgroup means (no CI)
)
writexl::write_xlsx(sheets, path = OUT_XLS)
message("[DONE] Fairness workbook written → ", OUT_XLS)

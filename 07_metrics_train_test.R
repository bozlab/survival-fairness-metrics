# ============================================================
# OS Performance — tau = ∞ (overall) AND fixed horizons {24,36,60} months
#   - Event/Time: Survival.months, event_status
#   - Predictors: 18 vars (NO Income_group, NO Urban_Rural)
#   - Selection: Univariable term p < .05 (fallback: TOP_K by smallest min-term-p)
#   - Metrics (single value per horizon): Harrell's C-index, IBS@tau, iAUC@tau
#   - BlackBoost: Breslow baseline adapter for predictSurvProb/predictRisk
#   - Output: one Excel with 5 sheets (tauInf, tau24, tau36, tau60, Combined)
# ============================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(broom)
  library(Hmisc)           # rcorr.cens (Harrell's C)
  library(riskRegression)  # Score(): Brier & AUC curves; predictRisk()
  library(pec)             # predictSurvProb()
  library(mboost)          # BlackBoost (CoxPH)
  library(randomForestSRC) # RSF
  library(writexl)
})

# ---------------------------
# Settings
# ---------------------------
set.seed(123)
ALPHA    <- 0.05   # univariable term p threshold
TOP_K    <- 5      # fallback when none pass ALPHA
GRID_N   <- 20     # # of time points for integration grid
BOOT_N   <- 1000    # bootstrap reps for CIs (↑ to 1000 for final)
HORIZONS <- c(24, 36, 60) # fixed horizons (months)

# ---------------------------
# Load data (OS)
# ---------------------------
train_data <- read.csv("train_data_seer.csv", check.names = TRUE)
test_data  <- read.csv("test_data_seer.csv",  check.names = TRUE)

# ---------------------------
# Guards: event/time types/coding
# ---------------------------
to01 <- function(x){
  if (is.factor(x)) x <- as.character(x)
  x <- as.integer(x); x[is.na(x)] <- 0L; x[x != 0L] <- 1L; x
}
if (!"event_status" %in% names(train_data))
  train_data$event_status <- ifelse(train_data$Vital.status.recode == "Dead", 1, 0)
if (!"event_status" %in% names(test_data))
  test_data$event_status  <- ifelse(test_data$Vital.status.recode == "Dead", 1, 0)

train_data$event_status    <- to01(train_data$event_status)
test_data$event_status     <- to01(test_data$event_status)
train_data$Survival.months <- as.numeric(train_data$Survival.months)
test_data$Survival.months  <- as.numeric(test_data$Survival.months)

# ---------------------------
# Predictors (18 vars; NO Income_group / NO Urban_Rural)
# ---------------------------
predictors <- c("Age_grouped","Race_grouped","Gender","Marital_grouped",
                "Primary_site_grouped","Histology_grouped","Grade_grouped","Laterality_grouped",
                "T_stage_recoded","N_stage_combined","M_stage_combined",
                "Surgery","Chemotherapy","Radiation",
                "Bone_met","Liver_met","Lung_met","Brain_met")

# Ensure factors (robust) and drop unused levels
for (v in predictors) {
  if (!is.null(train_data[[v]]) && !is.numeric(train_data[[v]])) train_data[[v]] <- factor(train_data[[v]])
  if (!is.null(test_data[[v]])  && !is.numeric(test_data[[v]]))  test_data[[v]]  <- factor(test_data[[v]])
}
train_data <- droplevels(train_data)
test_data  <- droplevels(test_data)

# ---------------------------
# Helper: align factor levels across splits
# ---------------------------
align_factors <- function(train, test, vars){
  for (v in vars) {
    if (is.character(train[[v]])) train[[v]] <- factor(train[[v]])
    if (is.character(test[[v]]))  test[[v]]  <- factor(test[[v]])
    if (is.factor(train[[v]]) || is.factor(test[[v]])) {
      lv <- union(levels(factor(train[[v]])), levels(factor(test[[v]])))
      train[[v]] <- factor(train[[v]], levels = lv)
      test[[v]]  <- factor(test[[v]],  levels = lv)
    }
  }
  list(train=train, test=test)
}

# ---------------------------
# Variable selection: univariable term p-values on TRAIN
# ---------------------------
surv_obj_tr <- with(train_data, Surv(Survival.months, event_status))
uni_results <- lapply(predictors, function(var) {
  fit <- coxph(as.formula(paste("surv_obj_tr ~", sprintf("`%s`", var))), data = train_data)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(variable = var)
})
uni_df <- bind_rows(uni_results)

# pick variables if any term p < ALPHA; otherwise TOP_K by smallest min-term-p
var_min_p <- uni_df %>% group_by(variable) %>%
  summarise(min_p = suppressWarnings(min(p.value, na.rm = TRUE)), .groups = "drop")
sig_vars <- var_min_p %>% filter(is.finite(min_p) & min_p < ALPHA) %>% pull(variable)
if (length(sig_vars) == 0) {
  message(sprintf("No variable passed p < %.2f. Falling back to TOP_%d by smallest p.", ALPHA, TOP_K))
  sig_vars <- var_min_p %>% arrange(min_p) %>% slice_head(n = min(TOP_K, n())) %>% pull(variable)
}

# Align factors for selected variables
al <- align_factors(train_data, test_data, sig_vars)
train_data <- al$train; test_data <- al$test

# Multivariable Cox (OS)
cox_formula <- as.formula(paste("Surv(Survival.months, event_status) ~", paste(sprintf("`%s`", sig_vars), collapse = " + ")))
print(cox_formula)
cox_model <- coxph(cox_formula, data = train_data, x = TRUE, y = TRUE, model = TRUE)

# ---------------------------
# RSF & BlackBoost (same covariates)
# ---------------------------
rsf_model <- randomForestSRC::rfsrc(cox_formula, data = train_data,
                                    ntree = 500, importance = "none", na.action = "na.impute")

bb_model  <- mboost::blackboost(cox_formula, data = train_data, family = mboost::CoxPH())
set.seed(123)
cvb <- cvrisk(bb_model, folds = cv(model.weights(bb_model), type = "kfold", B = 5))
mstop(bb_model) <- mstop(cvb)

# ---------------------------
# BlackBoost survival adapters (Breslow baseline + predict* generics)
# ---------------------------
.compute_bb_baseline <- function(model, time, status, newdata_for_eta){
  eta <- as.numeric(predict(model, newdata = newdata_for_eta, type = "link"))
  ord <- order(time, na.last = NA)
  time <- time[ord]; status <- status[ord]; eta <- eta[ord]
  ev_times <- sort(unique(time[status == 1]))
  H0 <- numeric(length(ev_times)); cum <- 0
  for (k in seq_along(ev_times)) {
    tk   <- ev_times[k]
    d_k  <- sum(status == 1 & time == tk)
    rset <- sum(exp(eta[time >= tk]))
    if (rset > 0) cum <- cum + d_k / rset
    H0[k] <- cum
  }
  list(times = ev_times, H0 = H0)
}
bb_baseline <- .compute_bb_baseline(bb_model,
                                    train_data$Survival.months,
                                    train_data$event_status,
                                    train_data)
attr(bb_model, "bb_baseline") <- bb_baseline

predictSurvProb.blackboost <- function(object, newdata, times, ...) {
  bl <- attr(object, "bb_baseline")
  if (is.null(bl) || is.null(bl$times) || is.null(bl$H0))
    stop("BlackBoost baseline not found. Attach 'bb_baseline'.")
  idx  <- findInterval(times, bl$times, left.open = FALSE)
  H0_t <- ifelse(idx == 0, 0, bl$H0[pmax(1, idx)])
  lp   <- as.numeric(predict(object, newdata = newdata, type = "link"))
  exp_eta <- exp(lp)
  N <- nrow(newdata); TT <- length(times)
  S <- matrix(NA_real_, nrow = N, ncol = TT)
  for (j in seq_len(TT)) S[, j] <- exp(- H0_t[j] * exp_eta)
  S
}
predictSurvProb.mboost <- function(object, newdata, times, ...) {
  predictSurvProb.blackboost(object, newdata = newdata, times = times, ...)
}
predictRisk.blackboost <- function(object, newdata, times, ...) {
  1 - predictSurvProb.blackboost(object, newdata = newdata, times = times, ...)
}
predictRisk.mboost <- function(object, newdata, times, ...) {
  1 - predictSurvProb.mboost(object, newdata = newdata, times = times, ...)
}

# ---------------------------
# Tau helpers + time grid + integration
# ---------------------------
# tau_inf := entire observed follow-up (min of max follow-up and max event time)
choose_tau_inf <- function(df){
  tmax <- suppressWarnings(max(df$Survival.months, na.rm = TRUE))
  tevt <- suppressWarnings(max(df$Survival.months[df$event_status == 1], na.rm = TRUE))
  tau  <- min(tmax, tevt)
  if (!is.finite(tau) || tau <= 0) tau <- 36
  tau
}
build_time_grid <- function(df, tau, n = GRID_N){
  t0 <- max(0.1, tau / (n * 2))     # avoid t=0
  seq(from = t0, to = tau, length.out = n)
}
trap_mean <- function(times, values){
  ok <- is.finite(times) & is.finite(values)
  times <- times[ok]; values <- values[ok]
  if (length(times) < 2) return(NA_real_)
  o <- order(times); times <- times[o]; values <- values[o]
  area <- sum(diff(times) * (head(values, -1) + tail(values, -1)) / 2)
  area / (max(times) - min(times))
}

curves_brier_auc <- function(data, model, label, times){
  sc <- tryCatch(
    riskRegression::Score(object = setNames(list(model), label),
                          formula = Surv(Survival.months, event_status) ~ 1,
                          data    = data,
                          metrics = c("Brier","AUC"),
                          times   = times,
                          cens.model   = "kaplan-meier",
                          split.method = "none",
                          conf.int     = FALSE),
    error = function(e) NULL
  )
  if (is.null(sc)) return(list(times = times, brier = rep(NA_real_, length(times)), auc = rep(NA_real_, length(times))))
  bd <- tryCatch(as.data.frame(sc$Brier$score), error = function(e) NULL)
  ad <- tryCatch(as.data.frame(sc$AUC$score),   error = function(e) NULL)
  bvec <- sapply(times, function(t) {
    if (is.null(bd)) return(NA_real_)
    row <- bd[bd$model == label & abs(bd$times - t) < 1e-8, , drop = FALSE]
    if (nrow(row)) as.numeric(row$Brier[1]) else NA_real_
  })
  avec <- sapply(times, function(t) {
    if (is.null(ad)) return(NA_real_)
    row <- ad[ad$model == label & abs(ad$times - t) < 1e-8, , drop = FALSE]
    if (nrow(row)) as.numeric(row$AUC[1]) else NA_real_
  })
  list(times = times, brier = bvec, auc = avec)
}

# single-run metrics on a dataset/model for a given tau
scalar_risk <- function(model, newdata, tau){
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "risk")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "lp")),   error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "link")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  # last resort: calibrated risk at tau
  mrk <- tryCatch({
    pr <- riskRegression::predictRisk(model, newdata = newdata, times = tau)
    as.numeric(pr[,1,drop=TRUE])
  }, error = function(e) rep(NA_real_, nrow(newdata)))
  mrk
}
metrics_once <- function(df, model, label, tau){
  r <- scalar_risk(model, df, tau)
  C <- tryCatch({
    as.numeric(Hmisc::rcorr.cens(-r, survival::Surv(df$Survival.months, df$event_status))["C Index"])
  }, error = function(e) NA_real_)
  times <- build_time_grid(df, tau, n = GRID_N)
  curv  <- curves_brier_auc(df, model, label, times)
  IBS   <- tryCatch(trap_mean(curv$times, curv$brier), error = function(e) NA_real_)
  iAUC  <- tryCatch(trap_mean(curv$times, curv$auc),   error = function(e) NA_real_)
  c(Cindex = C, IBS = IBS, iAUC = iAUC)
}

# bootstrap CIs
boot_metrics <- function(df, model, label, tau, n_boot = BOOT_N){
  mat <- matrix(NA_real_, nrow = n_boot, ncol = 3); colnames(mat) <- c("Cindex","IBS","iAUC")
  n <- nrow(df)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, replace = TRUE)
    boot <- df[idx, , drop = FALSE]
    if (nrow(boot) >= 30 && length(unique(boot$event_status)) >= 2) {
      mat[b, ] <- metrics_once(boot, model, label, tau)
    }
  }
  summarize <- function(v){
    v <- v[is.finite(v)]
    if (!length(v)) return(c(mean=NA, lower=NA, upper=NA, n=0))
    c(mean = mean(v), lower = as.numeric(quantile(v, 0.025)),
      upper = as.numeric(quantile(v, 0.975)), n = length(v))
  }
  list(Cindex = summarize(mat[, "Cindex"]),
       IBS    = summarize(mat[, "IBS"]),
       iAUC   = summarize(mat[, "iAUC"]))
}
fmt_ci3 <- function(m, lo, hi){
  if (any(!is.finite(c(m, lo, hi)))) return("NA [NA, NA]")
  sprintf("%.3f [%.3f, %.3f]", m, lo, hi)
}
row_for <- function(method, dataset, res, tag){
  data.frame(
    Horizon = tag,
    Method  = method,
    Dataset = dataset,
    `C-index (95% CI)` = fmt_ci3(res$Cindex["mean"], res$Cindex["lower"], res$Cindex["upper"]),
    `IBS (95% CI)`     = fmt_ci3(res$IBS["mean"],    res$IBS["lower"],    res$IBS["upper"]),
    `iAUC (95% CI)`    = fmt_ci3(res$iAUC["mean"],   res$iAUC["lower"],   res$iAUC["upper"]),
    check.names = FALSE
  )
}

# ---------------------------
# Evaluate tau = ∞ and fixed horizons (24, 36, 60)
# ---------------------------
# tau = ∞ (overall observed, computed per split)
choose_tau_inf <- function(df){
  tmax <- suppressWarnings(max(df$Survival.months, na.rm = TRUE))
  tevt <- suppressWarnings(max(df$Survival.months[df$event_status == 1], na.rm = TRUE))
  tau  <- min(tmax, tevt)
  if (!is.finite(tau) || tau <= 0) tau <- 36
  tau
}

tau_train_inf <- choose_tau_inf(train_data)
tau_test_inf  <- choose_tau_inf(test_data)
cat(sprintf("tau∞ — Train: %.2f months | Test: %.2f months\n", tau_train_inf, tau_test_inf))

cox_train_inf <- boot_metrics(train_data, cox_model, "Cox",        tau_train_inf, n_boot = BOOT_N)
cox_test_inf  <- boot_metrics(test_data,  cox_model, "Cox",        tau_test_inf,  n_boot = BOOT_N)
rsf_train_inf <- boot_metrics(train_data, rsf_model, "RSF",        tau_train_inf, n_boot = BOOT_N)
rsf_test_inf  <- boot_metrics(test_data,  rsf_model, "RSF",        tau_test_inf,  n_boot = BOOT_N)
bb_train_inf  <- boot_metrics(train_data, bb_model,  "BlackBoost", tau_train_inf, n_boot = BOOT_N)
bb_test_inf   <- boot_metrics(test_data,  bb_model,  "BlackBoost", tau_test_inf,  n_boot = BOOT_N)

table_inf <- bind_rows(
  row_for("Cox",        "Train", cox_train_inf, "tau=∞"),
  row_for("Cox",        "Test",  cox_test_inf,  "tau=∞"),
  row_for("RSF",        "Train", rsf_train_inf, "tau=∞"),
  row_for("RSF",        "Test",  rsf_test_inf,  "tau=∞"),
  row_for("BlackBoost", "Train", bb_train_inf,  "tau=∞"),
  row_for("BlackBoost", "Test",  bb_test_inf,   "tau=∞")
)

# Helper to evaluate a fixed horizon on both splits (with clipping to observed)
eval_fixed_tau <- function(tau_fixed){
  tau_tr <- min(choose_tau_inf(train_data), tau_fixed)
  tau_te <- min(choose_tau_inf(test_data),  tau_fixed)
  cat(sprintf("tau=%d — Train: %.2f | Test: %.2f\n", tau_fixed, tau_tr, tau_te))
  list(
    CoxTrain = boot_metrics(train_data, cox_model, "Cox",        tau_tr, n_boot = BOOT_N),
    CoxTest  = boot_metrics(test_data,  cox_model, "Cox",        tau_te, n_boot = BOOT_N),
    RSFTrain = boot_metrics(train_data, rsf_model, "RSF",        tau_tr, n_boot = BOOT_N),
    RSFTest  = boot_metrics(test_data,  rsf_model, "RSF",        tau_te, n_boot = BOOT_N),
    BBTrain  = boot_metrics(train_data, bb_model,  "BlackBoost", tau_tr, n_boot = BOOT_N),
    BBTest   = boot_metrics(test_data,  bb_model,  "BlackBoost", tau_te, n_boot = BOOT_N)
  )
}

# tau=24
res24 <- eval_fixed_tau(24)
table_24 <- bind_rows(
  row_for("Cox",        "Train", res24$CoxTrain, "tau=24"),
  row_for("Cox",        "Test",  res24$CoxTest,  "tau=24"),
  row_for("RSF",        "Train", res24$RSFTrain, "tau=24"),
  row_for("RSF",        "Test",  res24$RSFTest,  "tau=24"),
  row_for("BlackBoost", "Train", res24$BBTrain,  "tau=24"),
  row_for("BlackBoost", "Test",  res24$BBTest,   "tau=24")
)

# tau=36
res36 <- eval_fixed_tau(36)
table_36 <- bind_rows(
  row_for("Cox",        "Train", res36$CoxTrain, "tau=36"),
  row_for("Cox",        "Test",  res36$CoxTest,  "tau=36"),
  row_for("RSF",        "Train", res36$RSFTrain, "tau=36"),
  row_for("RSF",        "Test",  res36$RSFTest,  "tau=36"),
  row_for("BlackBoost", "Train", res36$BBTrain,  "tau=36"),
  row_for("BlackBoost", "Test",  res36$BBTest,   "tau=36")
)

# tau=60
res60 <- eval_fixed_tau(60)
table_60 <- bind_rows(
  row_for("Cox",        "Train", res60$CoxTrain, "tau=60"),
  row_for("Cox",        "Test",  res60$CoxTest,  "tau=60"),
  row_for("RSF",        "Train", res60$RSFTrain, "tau=60"),
  row_for("RSF",        "Test",  res60$RSFTest,  "tau=60"),
  row_for("BlackBoost", "Train", res60$BBTrain,  "tau=60"),
  row_for("BlackBoost", "Test",  res60$BBTest,   "tau=60")
)

# ---------------------------
# Save to Excel (5 sheets)
# ---------------------------
out_path <- "Table3_OS_Performance_tauInf_24_36_60.xlsx"
write_xlsx(
  list(
    tauInf   = table_inf   %>% select(-Horizon),
    tau24    = table_24    %>% select(-Horizon),
    tau36    = table_36    %>% select(-Horizon),
    tau60    = table_60    %>% select(-Horizon),
    Combined = bind_rows(table_inf, table_24, table_36, table_60)
  ),
  path = out_path
)
cat("Saved → ", out_path, "\n", sep = "")

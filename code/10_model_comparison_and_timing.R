suppressPackageStartupMessages({
  library(survival); library(dplyr); library(broom);  library(Hmisc)
  library(riskRegression); library(pec); library(mboost)
  library(randomForestSRC); library(writexl); library(parallel)
})

# --------------------------- Session design (macOS-friendly)
options(stringsAsFactors = FALSE, warn = 1)                     # flush warnings as they happen
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")        # avoid BLAS/OpenMP oversubscription
set.seed(123); RNGkind("L'Ecuyer-CMRG")
if (.Platform$OS.type == "unix") parallel::mc.reset.stream()   # deterministic mclapply streams on macOS
NCORES <- max(1, detectCores() - 1)
options(mc.cores = NCORES)

# Simple logging helpers
.ts   <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_i <- function(...) cat(sprintf("[%s] ", .ts()), sprintf(...), "\n", sep = "")

# --------------------------- NEW: Timing & size helpers
TIMES <- list()
time_it <- function(label, expr) {
  gc(); t0 <- proc.time()
  val <- force(expr)
  el <- as.numeric((proc.time() - t0)[["elapsed"]])
  TIMES[[label]] <<- el
  log_i("Timing: %s = %.3fs", label, el)
  return(val)
}
mb <- function(obj) as.numeric(utils::object.size(obj)) / (1024^2)  # MB
per_case_ms <- function(seconds, n) ifelse(n > 0, 1000 * seconds / n, NA_real_)

# --------------------------- Settings
ALPHA       <- 0.05
TOP_K       <- 5
GRID_N      <- 20          
BOOT_N      <- 1000        
FAIR_BOOT_N <- 1000       
HORIZONS    <- c(24,36,60)

# Fairness variables (not in model, only for evaluation)
FAIR_VARS <- c("Gender","Race_grouped","Age_grouped","Marital_grouped",
               "Income_group","Urban_Rural")
GAMMA_F   <- 0.01

# Intersectional: only 2-way
INTERSECT_MAX_K <- 2
MIN_GROUP_N     <- 25

# --------------------------- Paths
OUT_XLS_PERF <- "Table_OS_Performance.xlsx"   # Excel #1
OUT_XLS_FAIR <- "Table_OS_Fairness.xlsx"      # Excel #2
MODELS_RDS   <- "fitted_models_OS.rds"
TEST_RDS     <- "test_data_seer_test.rds"
TRAIN_CSV    <- "train_data_seer.csv"
TEST_CSV     <- "test_data_seer.csv"

# --------------------------- Load data
stopifnot(file.exists(TRAIN_CSV), file.exists(TEST_CSV))
train_data <- read.csv(TRAIN_CSV, check.names = TRUE)
test_data  <- read.csv(TEST_CSV,  check.names = TRUE)
log_i("Loaded data: train=%d rows, test=%d rows.", nrow(train_data), nrow(test_data))

# 🔎 Optional: downsample for smoke test (commented out for final)
# set.seed(42)
# train_data <- train_data[sample(nrow(train_data), 500), ]
# test_data  <- test_data[sample(nrow(test_data), 500), ]
# log_i("Downsampled for smoke test: train=%d, test=%d.", nrow(train_data), nrow(test_data))

# --------------------------- Outcome guards
to01 <- function(x){
  if (is.factor(x)) x <- as.character(x)
  x <- as.integer(x); x[is.na(x)] <- 0L; x[x != 0L] <- 1L; x
}
if (!"event_status" %in% names(train_data))
  train_data$event_status <- ifelse(train_data$Vital.status.recode=="Dead",1,0)
if (!"event_status" %in% names(test_data))
  test_data$event_status  <- ifelse(test_data$Vital.status.recode=="Dead",1,0)

train_data$event_status    <- to01(train_data$event_status)
test_data$event_status     <- to01(test_data$event_status)
train_data$Survival.months <- as.numeric(train_data$Survival.months)
test_data$Survival.months  <- as.numeric(test_data$Survival.months)

# --------------------------- Predictors (EXCLUDE Income_group & Urban_Rural from modeling)
predictors <- c(
  "Age_grouped","Race_grouped","Gender","Marital_grouped",
  "Primary_site_grouped","Histology_grouped","Grade_grouped","Laterality_grouped",
  "T_stage_recoded","N_stage_combined","M_stage_combined",
  "Surgery","Chemotherapy","Radiation",
  "Bone_met","Liver_met","Lung_met","Brain_met"
)

# Coerce to factors where appropriate and drop unused levels
for (v in predictors){
  if (!is.null(train_data[[v]]) && !is.numeric(train_data[[v]])) train_data[[v]] <- factor(train_data[[v]])
  if (!is.null(test_data[[v]])  && !is.numeric(test_data[[v]]))  test_data[[v]]  <- factor(test_data[[v]])
}
train_data <- droplevels(train_data); test_data <- droplevels(test_data)

# --------------------------- Factor alignment helper
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

# --------------------------- Response guard
complete_surv <- function(df){
  idx <- is.finite(df$Survival.months) & !is.na(df$event_status)
  df[idx,,drop=FALSE]
}

# --------------------------- Univariate screen -> sig_vars
surv_obj_tr <- with(train_data, Surv(Survival.months, event_status))
uni_results <- lapply(predictors, function(var){
  fit <- coxph(as.formula(paste("surv_obj_tr ~", sprintf("`%s`", var))), data = train_data)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(variable = var)
})
uni_df <- bind_rows(uni_results)
var_min_p <- uni_df %>%
  group_by(variable) %>%
  summarise(min_p = suppressWarnings(min(p.value, na.rm = TRUE)), .groups = "drop")

sig_vars <- var_min_p %>% filter(is.finite(min_p) & min_p < ALPHA) %>% pull(variable)
if (!length(sig_vars)){
  log_i("No variable passed p < %.2f. Falling back to TOP_%d by smallest p.", ALPHA, TOP_K)
  sig_vars <- var_min_p %>% arrange(min_p) %>% slice_head(n = min(TOP_K, n())) %>% pull(variable)
}
al <- align_factors(train_data, test_data, sig_vars); train_data <- al$train; test_data <- al$test

cox_formula <- as.formula(paste(
  "Surv(Survival.months, event_status) ~",
  paste(sprintf("`%s`", sig_vars), collapse = " + ")
))
log_i("Cox formula: %s", deparse(cox_formula))

# --------------------------- Models (with timing)
p <- length(sig_vars)

cox_model <- time_it("fit_cox", {
  coxph(cox_formula, data = train_data, x = TRUE, y = TRUE, model = TRUE)
})

rsf_model <- time_it("fit_rsf", {
  randomForestSRC::rfsrc(
    cox_formula, data = train_data,
    ntree = 1000, splitrule = "logrank",
    mtry = max(1, floor(sqrt(p))), nodesize = 15, nsplit = 0,
    importance = "none", na.action = "na.impute"
  )
})

bb_model <- time_it("fit_blackboost", {
  mboost::blackboost(cox_formula, data = train_data, family = mboost::CoxPH())
})
set.seed(123)
cvb <- time_it("cv_blackboost", { cvrisk(bb_model, folds = cv(model.weights(bb_model), type = "kfold", B = 5)) })
mstop(bb_model) <- mstop(cvb)
log_i("Models fitted.")

# --------------------------- Model sizes (MB)
MODEL_SIZES_MB <- c(
  Cox        = mb(cox_model),
  RSF        = mb(rsf_model),
  BlackBoost = mb(bb_model)
)
log_i("Model sizes (MB) — Cox=%.2f, RSF=%.2f, BlackBoost=%.2f",
      MODEL_SIZES_MB["Cox"], MODEL_SIZES_MB["RSF"], MODEL_SIZES_MB["BlackBoost"])

# --------------------------- BlackBoost baseline adapters (unchanged)
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
bb_baseline <- .compute_bb_baseline(bb_model, train_data$Survival.months, train_data$event_status, train_data)
attr(bb_model, "bb_baseline") <- bb_baseline

predictSurvProb.blackboost <- function(object, newdata, times, ...){
  bl <- attr(object, "bb_baseline")
  if (is.null(bl) || is.null(bl$times) || is.null(bl$H0)) stop("BlackBoost baseline not found. Attach 'bb_baseline'.")
  idx  <- findInterval(times, bl$times, left.open = FALSE)
  H0_t <- ifelse(idx == 0, 0, bl$H0[pmax(1, idx)])
  lp   <- as.numeric(predict(object, newdata = newdata, type = "link"))
  exp_eta <- exp(lp)
  N <- nrow(newdata); TT <- length(times); S <- matrix(NA_real_, nrow = N, ncol = TT)
  for (j in seq_len(TT)) S[, j] <- exp(- H0_t[j] * exp_eta)
  S
}
predictSurvProb.mboost <- function(object, newdata, times, ...) predictSurvProb.blackboost(object, newdata, times, ...)
predictRisk.blackboost <- function(object, newdata, times, ...) 1 - predictSurvProb.blackboost(object, newdata, times, ...)
predictRisk.mboost     <- function(object, newdata, times, ...) 1 - predictSurvProb.mboost(object, newdata, times, ...)

# --------------------------- RSF row-order-safe predictions (unchanged)
rsf_predict_surv_safe <- function(object, newdata, times){
  pr <- predict(object, newdata = newdata, na.action = "na.impute")
  tt <- pr$time.interest
  idx <- pmax(0L, findInterval(times, tt))
  S <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  for (j in seq_along(times)) if (idx[j] > 0L) S[, j] <- pr$survival[, idx[j]]
  S
}

# --------------------------- Unified survival matrix (unchanged)
surv_mat <- function(model, newdata, times){
  if (inherits(model, "coxph"))  return(pec::predictSurvProb(model, newdata = newdata, times = times))
  if (inherits(model, "mboost")) return(predictSurvProb(model, newdata = newdata, times = times))
  if (inherits(model, "rfsrc"))  return(rsf_predict_surv_safe(model, newdata = newdata, times = times))
  stop("Unsupported model class in surv_mat()")
}

# --------------------------- Tau & integration utils (unchanged)
choose_tau_inf <- function(df){
  tmax <- suppressWarnings(max(df$Survival.months, na.rm = TRUE))
  tevt <- suppressWarnings(max(df$Survival.months[df$event_status == 1], na.rm = TRUE))
  tau  <- min(tmax, tevt); if (!is.finite(tau) || tau <= 0) tau <- 36; tau
}
build_time_grid <- function(df, tau, n = GRID_N){
  t0 <- max(0.1, tau / (n * 2)); seq(from = t0, to = tau, length.out = n)
}
trap_mean <- function(times, values){
  ok <- is.finite(times) & is.finite(values)
  times <- times[ok]; values <- values[ok]
  if (length(times) < 2) return(NA_real_)
  o <- order(times); times <- times[o]; values <- values[o]
  area <- sum(diff(times) * (head(values, -1) + tail(values, -1)) / 2)
  area / (max(times) - min(times))
}

# --------------------------- Safe Score() wrapper (unchanged)
safe_Score <- function(risk_mat, data, times){
  ok <- rowSums(is.finite(risk_mat)) == ncol(risk_mat)
  if (!all(ok)){
    data     <- data[ok,,drop=FALSE]
    risk_mat <- risk_mat[ok,,drop=FALSE]
  }
  if (nrow(data) < 30 || length(unique(data$event_status)) < 2)
    return(list(Brier = rep(NA_real_, length(times)), AUC = rep(NA_real_, length(times))))
  res <- try({
    sc <- riskRegression::Score(
      object       = list(label = risk_mat),
      formula      = Surv(Survival.months, event_status) ~ 1,
      data         = data,
      metrics      = c("Brier","AUC"),
      times        = times,
      cens.model   = "kaplan-meier",
      split.method = "none",
      conf.int     = FALSE
    )
    bd <- as.data.frame(sc$Brier$score)
    ad <- as.data.frame(sc$AUC$score)
    B  <- sapply(times, function(t) as.numeric(bd$Brier[bd$model=="label" & abs(bd$times - t) < 1e-8][1]))
    A  <- sapply(times, function(t) as.numeric(ad$AUC[  ad$model=="label" & abs(ad$times - t) < 1e-8][1]))
    list(Brier = B, AUC = A)
  }, silent = TRUE)
  if (inherits(res, "try-error")) {
    log_i("Score() returned error; filled NA for Brier/AUC.")
    return(list(Brier = rep(NA_real_, length(times)), AUC = rep(NA_real_, length(times))))
  }
  res
}

# --------------------------- Curves (Brier/AUC) — robust design (unchanged)
curves_brier_auc <- function(data, model, label, times){
  df <- complete_surv(data)
  if (nrow(df) < 30 || length(unique(df$event_status)) < 2)
    return(list(times = times, brier = rep(NA_real_, length(times)), auc = rep(NA_real_, length(times))))
  maxfup <- suppressWarnings(max(df$Survival.months, na.rm = TRUE))
  times  <- sort(unique(pmin(times[times > 0 & is.finite(times)], maxfup)))
  S <- surv_mat(model, df, times)
  risk_mat <- 1 - S
  colnames(risk_mat) <- as.character(times)
  sc <- safe_Score(risk_mat, df, times)
  list(times = times, brier = sc$Brier, auc = sc$AUC)
}

# --------------------------- Scalar risk for C-index (unchanged)
scalar_risk <- function(model, newdata, tau){
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "risk")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "lp")),   error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  r <- tryCatch(as.numeric(predict(model, newdata = newdata, type = "link")), error = function(e) NULL)
  if (!is.null(r) && all(is.finite(r))) return(r)
  if (inherits(model,"rfsrc")){
    S <- rsf_predict_surv_safe(model, newdata, times = tau)
    return(1 - as.numeric(S[,1,drop=TRUE]))
  }
  mrk <- tryCatch({
    pr <- riskRegression::predictRisk(model, newdata = newdata, times = tau)
    as.numeric(pr[,1,drop=TRUE])
  }, error = function(e) rep(NA_real_, nrow(newdata)))
  mrk
}

# --------------------------- Single-run metrics (unchanged)
metrics_once <- function(df, model, label, tau){
  df <- complete_surv(df)
  if (nrow(df) < 30 || length(unique(df$event_status)) < 2)
    return(c(Cindex = NA_real_, IBS = NA_real_, iAUC = NA_real_))
  r <- scalar_risk(model, df, tau)
  okr <- is.finite(r)
  C <- if (sum(okr) < 30 || length(unique(df$event_status[okr])) < 2) NA_real_ else
    tryCatch(as.numeric(Hmisc::rcorr.cens(-r[okr], Surv(df$Survival.months[okr], df$event_status[okr]))["C Index"]),
             error = function(e) NA_real_)
  times <- build_time_grid(df, tau, n = GRID_N)
  curv  <- curves_brier_auc(df, model, label, times)
  IBS   <- tryCatch(trap_mean(curv$times, curv$brier), error = function(e) NA_real_)
  iAUC  <- tryCatch(trap_mean(curv$times, curv$auc),   error = function(e) NA_real_)
  c(Cindex = C, IBS = IBS, iAUC = iAUC)
}

# --------------------------- Pretty format (unchanged)
fmt_ci3 <- function(m, lo, hi){
  if (any(!is.finite(c(m, lo, hi)))) return("NA [NA, NA]")
  sprintf("%.3f [%.3f, %.3f]", m, lo, hi)
}
row_for <- function(method, dataset, res, tag){
  data.frame(
    Horizon = tag, Method = method, Dataset = dataset,
    `C-index (95% CI)` = fmt_ci3(res$Cindex["mean"], res$Cindex["lower"], res$Cindex["upper"]),
    `IBS (95% CI)`     = fmt_ci3(res$IBS["mean"],    res$IBS["lower"],    res$IBS["upper"]),
    `iAUC (95% CI)`    = fmt_ci3(res$iAUC["mean"],   res$iAUC["lower"],   res$iAUC["upper"]),
    check.names = FALSE
  )
}

# --------------------------- Bootstrap engine with progress (unchanged)
boot_metrics <- function(df, model, label, tau, n_boot = BOOT_N, ncores = NCORES){
  df <- df
  boot_core <- function(b){
    boot <- df[sample.int(nrow(df), replace = TRUE), , drop = FALSE]
    boot <- complete_surv(boot)
    if (nrow(boot) < 30 || length(unique(boot$event_status)) < 2)
      return(c(Cindex = NA_real_, IBS = NA_real_, iAUC = NA_real_))
    metrics_once(boot, model, label, tau)
  }
  pb <- txtProgressBar(min = 0, max = 100, style = 3); on.exit(close(pb), add = TRUE)
  step <- max(1L, floor(n_boot/100L)); ids <- split(seq_len(n_boot), ceiling(seq_len(n_boot)/step))
  res_list <- list(); k_done <- 0L
  for (chunk in ids){
    chunk_res <- mclapply(chunk, boot_core, mc.cores = ncores)
    res_list  <- c(res_list, chunk_res); k_done <- k_done + length(chunk)
    setTxtProgressBar(pb, floor(100*k_done/n_boot))
  }
  mat <- do.call(rbind, res_list); colnames(mat) <- c("Cindex","IBS","iAUC")
  summarize <- function(v){
    v <- v[is.finite(v)]; if (!length(v)) return(c(mean = NA, lower = NA, upper = NA, n = 0))
    c(mean = mean(v), lower = as.numeric(quantile(v, 0.025)),
      upper = as.numeric(quantile(v, 0.975)), n = length(v))
  }
  list(Cindex = summarize(mat[, "Cindex"]),
       IBS    = summarize(mat[, "IBS"]),
       iAUC   = summarize(mat[, "iAUC"]))
}

# ============================================================
# OVERALL evaluation (tau = ∞ and fixed horizons)
# ============================================================
tau_train_inf <- choose_tau_inf(train_data)
tau_test_inf  <- choose_tau_inf(test_data)
log_i("Tau∞ — Train: %.2f months | Test: %.2f months", tau_train_inf, tau_test_inf)

# ---- NEW: One-pass metrics timing on TEST (tau∞)
METR_COX <- time_it("metrics_once_cox_test",  metrics_once(test_data, cox_model, "Cox",        tau_test_inf))
METR_RSF <- time_it("metrics_once_rsf_test",  metrics_once(test_data, rsf_model, "RSF",        tau_test_inf))
METR_BB  <- time_it("metrics_once_bb_test",   metrics_once(test_data, bb_model,  "BlackBoost", tau_test_inf))

# ---- NEW: Inference timing (survival matrices on TEST for grid)
times_inf <- build_time_grid(test_data, tau_test_inf, n = GRID_N)
S_cox <- time_it("infer_cox_survmat",  surv_mat(cox_model, test_data, times_inf))
S_rsf <- time_it("infer_rsf_survmat",  surv_mat(rsf_model, test_data, times_inf))
S_bb  <- time_it("infer_bb_survmat",   surv_mat(bb_model,  test_data, times_inf))

INF_MS_PER_CASE <- c(
  Cox        = per_case_ms(TIMES[["infer_cox_survmat"]], nrow(test_data)),
  RSF        = per_case_ms(TIMES[["infer_rsf_survmat"]], nrow(test_data)),
  BlackBoost = per_case_ms(TIMES[["infer_bb_survmat"]],  nrow(test_data))
)
log_i("Per-case inference (ms) — Cox=%.2f, RSF=%.2f, BlackBoost=%.2f",
      INF_MS_PER_CASE["Cox"], INF_MS_PER_CASE["RSF"], INF_MS_PER_CASE["BlackBoost"])

# ---- Bootstrap (Train/Test; keeping your original structure)
cox_train_inf <- boot_metrics(train_data, cox_model, "Cox",        tau_train_inf)
cox_test_inf  <- time_it("boot_cox_test",  boot_metrics(test_data,  cox_model, "Cox",        tau_test_inf, n_boot = BOOT_N, ncores = NCORES))
rsf_train_inf <- boot_metrics(train_data, rsf_model, "RSF",        tau_train_inf)
rsf_test_inf  <- time_it("boot_rsf_test",  boot_metrics(test_data,  rsf_model, "RSF",        tau_test_inf, n_boot = BOOT_N, ncores = NCORES))
bb_train_inf  <- boot_metrics(train_data, bb_model,  "BlackBoost", tau_train_inf)
bb_test_inf   <- time_it("boot_bb_test",   boot_metrics(test_data,  bb_model,  "BlackBoost", tau_test_inf, n_boot = BOOT_N, ncores = NCORES))

BOOT_PER_IT_SEC <- c(
  Cox        = TIMES[["boot_cox_test"]] / BOOT_N,
  RSF        = TIMES[["boot_rsf_test"]] / BOOT_N,
  BlackBoost = TIMES[["boot_bb_test"]]  / BOOT_N
)
log_i("Bootstrap per-iter (s) — Cox=%.3f, RSF=%.3f, BlackBoost=%.3f",
      BOOT_PER_IT_SEC["Cox"], BOOT_PER_IT_SEC["RSF"], BOOT_PER_IT_SEC["BlackBoost"])

table_inf <- bind_rows(
  row_for("Cox","Train",cox_train_inf,"tau=∞"),
  row_for("Cox","Test", cox_test_inf,"tau=∞"),
  row_for("RSF","Train",rsf_train_inf,"tau=∞"),
  row_for("RSF","Test", rsf_test_inf,"tau=∞"),
  row_for("BlackBoost","Train",bb_train_inf,"tau=∞"),
  row_for("BlackBoost","Test", bb_test_inf,"tau=∞")
)

# Keep your fixed-tau evaluator (unchanged)
eval_fixed_tau <- function(tau_fixed){
  tau_tr <- min(choose_tau_inf(train_data), tau_fixed)
  tau_te <- min(choose_tau_inf(test_data),  tau_fixed)
  log_i("Tau=%d months — Train: %.2f | Test: %.2f", tau_fixed, tau_tr, tau_te)
  list(
    CoxTrain = boot_metrics(train_data, cox_model, "Cox",        tau_tr),
    CoxTest  = boot_metrics(test_data,  cox_model, "Cox",        tau_te),
    RSFTrain = boot_metrics(train_data, rsf_model, "RSF",        tau_tr),
    RSFTest  = boot_metrics(test_data,  rsf_model, "RSF",        tau_te),
    BBTrain  = boot_metrics(train_data, bb_model,  "BlackBoost", tau_tr),
    BBTest   = boot_metrics(test_data,  bb_model,  "BlackBoost", tau_te)
  )
}

# --------------------------- NEW: Publish-ready timing table
timing_table <- data.frame(
  Method = c("Cox", "RSF", "BlackBoost"),
  Fit_sec = c(TIMES[["fit_cox"]], TIMES[["fit_rsf"]], TIMES[["fit_blackboost"]] + TIMES[["cv_blackboost"]]),
  Inference_total_sec = c(TIMES[["infer_cox_survmat"]], TIMES[["infer_rsf_survmat"]], TIMES[["infer_bb_survmat"]]),
  Inference_per_case_ms = as.numeric(INF_MS_PER_CASE[c("Cox","RSF","BlackBoost")]),
  Onepass_metrics_sec = c(TIMES[["metrics_once_cox_test"]], TIMES[["metrics_once_rsf_test"]], TIMES[["metrics_once_bb_test"]]),
  Bootstrap_total_sec = c(TIMES[["boot_cox_test"]], TIMES[["boot_rsf_test"]], TIMES[["boot_bb_test"]]),
  Bootstrap_per_iter_sec = as.numeric(BOOT_PER_IT_SEC[c("Cox","RSF","BlackBoost")]),
  Model_size_MB = as.numeric(MODEL_SIZES_MB[c("Cox","RSF","BlackBoost")]),
  stringsAsFactors = FALSE
)
log_i("Timing table ready."); print(timing_table)

# ---- (Optional) Write to Excel alongside performance table
# writexl::write_xlsx(list(
#   "Performance (tau∞)" = table_inf,
#   "Timing (tau∞ Test)" = timing_table
# ), OUT_XLS_PERF)



# extract bootstrap means
boot_means <- data.frame(
  Method = c("Cox","RSF","BlackBoost"),
  Cindex_boot = c(cox_test_inf$Cindex["mean"],
                  rsf_test_inf$Cindex["mean"],
                  bb_test_inf$Cindex["mean"]),
  IBS_boot    = c(cox_test_inf$IBS["mean"],
                  rsf_test_inf$IBS["mean"],
                  bb_test_inf$IBS["mean"]),
  iAUC_boot   = c(cox_test_inf$iAUC["mean"],
                  rsf_test_inf$iAUC["mean"],
                  bb_test_inf$iAUC["mean"])
)

compare_table <- merge(onepass_print, boot_means, by="Method")
print(compare_table, row.names = FALSE)
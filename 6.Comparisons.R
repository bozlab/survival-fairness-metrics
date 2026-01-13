

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

# --------------------------- Settings
ALPHA       <- 0.05
TOP_K       <- 5
GRID_N      <- 20          # time grid for IBS/iAUC integration (20–30 recommended)
BOOT_N      <- 1000        # final: 1000 (use 2/50 for quick trials)
FAIR_BOOT_N <- 1000        # final: 1000
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

# --------------------------- Models
t0 <- proc.time()
cox_model <- coxph(cox_formula, data = train_data, x = TRUE, y = TRUE, model = TRUE)

p <- length(sig_vars)
rsf_model <- randomForestSRC::rfsrc(
  cox_formula, data = train_data,
  ntree = 1000, splitrule = "logrank",
  mtry = max(1, floor(sqrt(p))), nodesize = 15, nsplit = 0,
  importance = "none", na.action = "na.impute"
)

bb_model <- mboost::blackboost(cox_formula, data = train_data, family = mboost::CoxPH())
set.seed(123)
cvb <- cvrisk(bb_model, folds = cv(model.weights(bb_model), type = "kfold", B = 5))
mstop(bb_model) <- mstop(cvb)
log_i("Models fitted. Elapsed: %.1fs", (proc.time() - t0)[["elapsed"]])

# --------------------------- BlackBoost baseline adapters
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

# --------------------------- RSF row-order-safe predictions
rsf_predict_surv_safe <- function(object, newdata, times){
  pr <- predict(object, newdata = newdata, na.action = "na.impute")
  tt <- pr$time.interest
  idx <- pmax(0L, findInterval(times, tt))
  S <- matrix(1, nrow = nrow(newdata), ncol = length(times))
  for (j in seq_along(times)) if (idx[j] > 0L) S[, j] <- pr$survival[, idx[j]]
  S
}

# --------------------------- Unified survival matrix
surv_mat <- function(model, newdata, times){
  if (inherits(model, "coxph"))  return(pec::predictSurvProb(model, newdata = newdata, times = times))
  if (inherits(model, "mboost")) return(predictSurvProb(model, newdata = newdata, times = times))
  if (inherits(model, "rfsrc"))  return(rsf_predict_surv_safe(model, newdata = newdata, times = times))
  stop("Unsupported model class in surv_mat()")
}

# --------------------------- Tau & integration utils
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

# --------------------------- Safe Score() wrapper (design fix)
safe_Score <- function(risk_mat, data, times){
  # Drop rows with any NA in risk predictions to avoid "Missing values in predicted risk"
  ok <- rowSums(is.finite(risk_mat)) == ncol(risk_mat)
  if (!all(ok)){
    data     <- data[ok,,drop=FALSE]
    risk_mat <- risk_mat[ok,,drop=FALSE]
  }
  if (nrow(data) < 30 || length(unique(data$event_status)) < 2)
    return(list(Brier = rep(NA_real_, length(times)), AUC = rep(NA_real_, length(times))))
  # Try Score(); if it fails, return NA vectors (never stop the pipeline)
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

# --------------------------- Curves (Brier/AUC) — robust design
curves_brier_auc <- function(data, model, label, times){
  df <- complete_surv(data)
  if (nrow(df) < 30 || length(unique(df$event_status)) < 2)
    return(list(times = times, brier = rep(NA_real_, length(times)), auc = rep(NA_real_, length(times))))
  # Clamp times to observed follow-up
  maxfup <- suppressWarnings(max(df$Survival.months, na.rm = TRUE))
  times  <- sort(unique(pmin(times[times > 0 & is.finite(times)], maxfup)))
  # Predict survival -> risks
  S <- surv_mat(model, df, times)
  risk_mat <- 1 - S
  colnames(risk_mat) <- as.character(times)
  # Safe scoring
  sc <- safe_Score(risk_mat, df, times)
  list(times = times, brier = sc$Brier, auc = sc$AUC)
}

# --------------------------- Scalar risk for C-index (fallbacks)
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

# --------------------------- Single-run metrics (robust to NA risks)
metrics_once <- function(df, model, label, tau){
  df <- complete_surv(df)
  if (nrow(df) < 30 || length(unique(df$event_status)) < 2)
    return(c(Cindex = NA_real_, IBS = NA_real_, iAUC = NA_real_))
  # C-index: drop NA risks
  r <- scalar_risk(model, df, tau)
  okr <- is.finite(r)
  C <- if (sum(okr) < 30 || length(unique(df$event_status[okr])) < 2) NA_real_ else
    tryCatch(as.numeric(Hmisc::rcorr.cens(-r[okr], Surv(df$Survival.months[okr], df$event_status[okr]))["C Index"]),
             error = function(e) NA_real_)
  # IBS / iAUC via curves
  times <- build_time_grid(df, tau, n = GRID_N)
  curv  <- curves_brier_auc(df, model, label, times)
  IBS   <- tryCatch(trap_mean(curv$times, curv$brier), error = function(e) NA_real_)
  iAUC  <- tryCatch(trap_mean(curv$times, curv$auc),   error = function(e) NA_real_)
  c(Cindex = C, IBS = IBS, iAUC = iAUC)
}

# --------------------------- Pretty format
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

# --------------------------- Bootstrap engine with progress
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

cox_train_inf <- boot_metrics(train_data, cox_model, "Cox",        tau_train_inf)
cox_test_inf  <- boot_metrics(test_data,  cox_model, "Cox",        tau_test_inf)
rsf_train_inf <- boot_metrics(train_data, rsf_model, "RSF",        tau_train_inf)
rsf_test_inf  <- boot_metrics(test_data,  rsf_model, "RSF",        tau_test_inf)
bb_train_inf  <- boot_metrics(train_data, bb_model,  "BlackBoost", tau_train_inf)
bb_test_inf   <- boot_metrics(test_data,  bb_model,  "BlackBoost", tau_test_inf)

table_inf <- bind_rows(
  row_for("Cox","Train",cox_train_inf,"tau=∞"),
  row_for("Cox","Test", cox_test_inf,"tau=∞"),
  row_for("RSF","Train",rsf_train_inf,"tau=∞"),
  row_for("RSF","Test", rsf_test_inf,"tau=∞"),
  row_for("BlackBoost","Train",bb_train_inf,"tau=∞"),
  row_for("BlackBoost","Test", bb_test_inf,"tau=∞")
)

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
res24 <- eval_fixed_tau(24); res36 <- eval_fixed_tau(36); res60 <- eval_fixed_tau(60)

table_24 <- bind_rows(
  row_for("Cox","Train",res24$CoxTrain,"tau=24"), row_for("Cox","Test", res24$CoxTest,"tau=24"),
  row_for("RSF","Train",res24$RSFTrain,"tau=24"), row_for("RSF","Test", res24$RSFTest,"tau=24"),
  row_for("BlackBoost","Train",res24$BBTrain,"tau=24"), row_for("BlackBoost","Test", res24$BBTest,"tau=24")
)
table_36 <- bind_rows(
  row_for("Cox","Train",res36$CoxTrain,"tau=36"), row_for("Cox","Test", res36$CoxTest,"tau=36"),
  row_for("RSF","Train",res36$RSFTrain,"tau=36"), row_for("RSF","Test", res36$RSFTest,"tau=36"),
  row_for("BlackBoost","Train",res36$BBTrain,"tau=36"), row_for("BlackBoost","Test", res36$BBTest,"tau=36")
)
table_60 <- bind_rows(
  row_for("Cox","Train",res60$CoxTrain,"tau=60"), row_for("Cox","Test", res60$CoxTest,"tau=60"),
  row_for("RSF","Train",res60$RSFTrain,"tau=60"), row_for("RSF","Test", res60$RSFTest,"tau=60"),
  row_for("BlackBoost","Train",res60$BBTrain,"tau=60"), row_for("BlackBoost","Test", res60$BBTest,"tau=60")
)

# --------------------------- WRITE EXCEL #1 (Performance)
write_xlsx(
  list(
    Overall_tauInf   = table_inf   %>% dplyr::select(-Horizon),
    Overall_tau24    = table_24    %>% dplyr::select(-Horizon),
    Overall_tau36    = table_36    %>% dplyr::select(-Horizon),
    Overall_tau60    = table_60    %>% dplyr::select(-Horizon),
    Overall_Combined = dplyr::bind_rows(table_inf, table_24, table_36, table_60)
  ),
  path = OUT_XLS_PERF
)
log_i("✅ Saved Performance tables → %s", OUT_XLS_PERF)

# =========================== FAIRNESS BLOCK (TEST only) ===========================
as_factor_safe <- function(df, var){
  if (!var %in% names(df)) return(df)
  if (!is.numeric(df[[var]])) df[[var]] <- factor(df[[var]])
  df
}
boot_row <- function(model, name, df, tau, tag){
  res <- boot_metrics(df, model, name, tau, n_boot = BOOT_N, ncores = NCORES)
  data.frame(Horizon = tag, Method = name, Dataset = "Test",
             `C-index (95% CI)` = fmt_ci3(res$Cindex["mean"],res$Cindex["lower"],res$Cindex["upper"]),
             `IBS (95% CI)`     = fmt_ci3(res$IBS["mean"],   res$IBS["lower"],   res$IBS["upper"]),
             `iAUC (95% CI)`    = fmt_ci3(res$iAUC["mean"],  res$iAUC["lower"],  res$iAUC["upper"]),
             check.names = FALSE)
}
fairness_for_var <- function(var, df_test, models, taus){
  if (!var %in% names(df_test)){ log_i("[FAIRNESS] Skipping '%s' (not in data).", var); return(NULL) }
  df_test <- as_factor_safe(df_test, var)
  lv <- levels(factor(df_test[[var]]))
  out <- list()
  for (L in lv){
    sub <- df_test[df_test[[var]] == L & !is.na(df_test[[var]]), , drop = FALSE]
    sub <- complete_surv(sub)
    if (nrow(sub) < 50 || length(unique(sub$event_status)) < 2){
      log_i("[FAIRNESS] %s=%s skipped (n<50 or single event class).", var, L); next
    }
    for (tau in taus){
      tau_use <- if (is.infinite(tau)) choose_tau_inf(sub) else min(choose_tau_inf(sub), tau)
      tag <- if (is.infinite(tau)) "tau=∞" else paste0("tau=", tau)
      for (mname in names(models)){
        out[[length(out) + 1]] <- cbind(Group = var, Level = L,
                                        boot_row(models[[mname]], mname, sub, tau_use, tag))
      }
    }
  }
  if (length(out)) bind_rows(out) else NULL
}
models_list <- list(Cox = cox_model, RSF = rsf_model, BlackBoost = bb_model)
taus_all    <- c(Inf, HORIZONS)
fair_tables <- lapply(FAIR_VARS, function(g) fairness_for_var(g, test_data, models_list, taus_all))
fairness_combined <- bind_rows(Filter(Negate(is.null), fair_tables))

# ---- CI / FCI / FCG helpers
group_CF <- function(risk, time, status, grp){
  levs <- levels(factor(grp))
  sapply(levs, function(L){
    idx <- which(grp == L)
    if (length(idx) < 30 || length(unique(status[idx])) < 2) return(NA_real_)
    as.numeric(Hmisc::rcorr.cens(-risk[idx], Surv(time[idx], status[idx]))["C Index"])
  })
}
compute_CI_point <- function(model, df, sens_var, times_for_prob){
  df <- complete_surv(df)
  A  <- factor(df[[sens_var]]); Tm <- df$Survival.months; St <- df$event_status
  if (inherits(model,"coxph") || inherits(model,"mboost")){
    r  <- scalar_risk(model, df, tau = median(Tm, na.rm = TRUE))
    ok <- is.finite(r); if (!any(ok)) return(NA_real_)
    CF <- group_CF(r[ok], Tm[ok], St[ok], A[ok])
    return(100 * max(abs(outer(CF, CF, `-`)), na.rm = TRUE))
  } else {
    S <- surv_mat(model, df, times_for_prob)
    cis <- vapply(seq_along(times_for_prob), function(j){
      rj <- 1 - S[, j]; ok <- is.finite(rj)
      CFj <- group_CF(rj[ok], Tm[ok], St[ok], A[ok])
      max(abs(outer(CFj, CFj, `-`)), na.rm = TRUE)
    }, numeric(1))
    return(100 * mean(cis, na.rm = TRUE))
  }
}
design_matrix <- function(df){
  X <- model.matrix(cox_formula, data = df)
  X <- X[, colnames(X) != "(Intercept)", drop = TRUE]
  scale(X)
}
pair_stat_one_t <- function(S, time, status, X, gamma, sens = NULL){
  Nc  <- which(status == 0)
  Nuc <- which(status == 1)
  if (!length(Nc) || !length(Nuc)) return(list(FCI = NA, FCG = NA))
  acc_FCI <- 0; cnt_FCI <- 0; acc_FCG <- 0; cnt_FCG <- 0
  for (i in Nc){
    validJ <- Nuc[ time[Nuc] >= time[i] ]
    if (length(validJ)){
      Si <- S[i];  Xi <- X[i,,drop=FALSE]
      Sj <- S[validJ]; Xj <- X[validJ,,drop=FALSE]
      dX <- sqrt(rowSums((sweep(Xj, 2, Xi))^2))
      term <- pmax(abs(Si - Sj) - gamma * dX, 0)
      acc_FCI <- acc_FCI + sum(term);        cnt_FCI <- cnt_FCI + length(term)
      if (!is.null(sens)){
        same <- which(sens[validJ] == sens[i])
        if (length(same)){ acc_FCG <- acc_FCG + sum(term[same]); cnt_FCG <- cnt_FCG + length(same) }
      }
    }
  }
  FCI <- if (cnt_FCI > 0) acc_FCI / (length(Nc) * length(Nuc)) else NA
  FCG <- if (cnt_FCG > 0) acc_FCG / (length(Nc) * length(Nuc)) else NA
  list(FCI = FCI, FCG = FCG)
}
compute_FCI_FCG <- function(model, df, sens_var, times, gamma){
  df2 <- complete_surv(df)
  if (nrow(df2) < 30) return(c(FCI = NA, FCG = NA))
  S    <- surv_mat(model, df2, times)
  X    <- design_matrix(df2)
  sens <- factor(df2[[sens_var]])
  t_res <- lapply(seq_along(times), function(j) pair_stat_one_t(S[, j], df2$Survival.months, df2$event_status, X, gamma, sens))
  FCI_mean <- mean(sapply(t_res, `[[`, "FCI"), na.rm = TRUE)
  FCG_mean <- mean(sapply(t_res, `[[`, "FCG"), na.rm = TRUE)
  c(FCI = FCI_mean, FCG = FCG_mean)
}

# ---- Overall fairness (point)
times_prob <- as.numeric(quantile(complete_surv(test_data)$Survival.months, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
fairness_point_tbl <- {
  rows <- list()
  for (v in FAIR_VARS) if (v %in% names(test_data)) for (mname in names(models_list)) {
    ci_val <- compute_CI_point(models_list[[mname]], test_data, v, times_prob)
    ff     <- compute_FCI_FCG(models_list[[mname]], test_data, v, times_prob, GAMMA_F)
    rows[[length(rows) + 1]] <- data.frame(Sensitive = v, Method = mname,
                                           CI_percent = round(ci_val, 3),
                                           F_CI = round(ff["FCI"], 4),
                                           F_CG = round(ff["FCG"], 4),
                                           check.names = FALSE)
  }
  dplyr::bind_rows(rows)
}

# ---- Overall fairness (bootstrap CIs)
fairness_boot_tbl <- {
  pb <- txtProgressBar(min = 0, max = 100, style = 3); on.exit(close(pb), add = TRUE)
  step <- max(1L, floor(FAIR_BOOT_N/100L))
  ids  <- split(seq_len(FAIR_BOOT_N), ceiling(seq_len(FAIR_BOOT_N)/step))
  res_all <- list(); done <- 0L
  for (chunk in ids) {
    chunk_res <- mclapply(chunk, function(b){
      idx  <- sample.int(nrow(test_data), replace = TRUE)
      boot <- test_data[idx,,drop=FALSE]
      rows <- list()
      for (v in FAIR_VARS) if (v %in% names(boot)) for (mname in names(models_list)) {
        ci_val <- compute_CI_point(models_list[[mname]], boot, v, times_prob)
        ff     <- compute_FCI_FCG(models_list[[mname]], boot, v, times_prob, GAMMA_F)
        rows[[length(rows) + 1]] <- data.frame(Sensitive = v, Method = mname,
                                               CI_percent = ci_val, F_CI = ff["FCI"], F_CG = ff["FCG"])
      }
      dplyr::bind_rows(rows)
    }, mc.cores = NCORES)
    res_all <- c(res_all, chunk_res); done <- done + length(chunk)
    setTxtProgressBar(pb, floor(100*done/FAIR_BOOT_N))
  }
  big <- dplyr::bind_rows(Map(function(df, rep_id) cbind(df, .rep = rep_id), res_all, seq_along(res_all)))
  big %>%
    dplyr::group_by(Sensitive, Method) %>%
    dplyr::summarise(
      CI_mean   = mean(CI_percent, na.rm = TRUE),
      CI_lower  = quantile(CI_percent, 0.025, na.rm = TRUE),
      CI_upper  = quantile(CI_percent, 0.975, na.rm = TRUE),
      CI_se     = sd(CI_percent, na.rm = TRUE),
      CI_n      = sum(is.finite(CI_percent)),
      FCI_mean  = mean(F_CI, na.rm = TRUE),
      FCI_lower = quantile(F_CI, 0.025, na.rm = TRUE),
      FCI_upper = quantile(F_CI, 0.975, na.rm = TRUE),
      FCI_se    = sd(F_CI, na.rm = TRUE),
      FCI_n     = sum(is.finite(F_CI)),
      FCG_mean  = mean(F_CG, na.rm = TRUE),
      FCG_lower = quantile(F_CG, 0.025, na.rm = TRUE),
      FCG_upper = quantile(F_CG, 0.975, na.rm = TRUE),
      FCG_se    = sd(F_CG, na.rm = TRUE),
      FCG_n     = sum(is.finite(F_CG)),
      .groups="drop"
    ) %>%
    dplyr::mutate(
      `CI(%) [95%CI]` = sprintf("%.3f [%.3f, %.3f]", CI_mean, CI_lower, CI_upper),
      `F_CI [95%CI]`  = sprintf("%.4f [%.4f, %.4f]", FCI_mean, FCI_lower, FCI_upper),
      `F_CG [95%CI]`  = sprintf("%.4f [%.4f, %.4f]", FCG_mean, FCG_lower, FCG_upper)
    )
}

# ---- 2-way intersectional F_cap (point + bootstrap)
build_intersect_sets <- function(vars, kmax = 2) {
  sets <- list(); for (k in 2:kmax) sets <- c(sets, combn(vars, k, simplify = FALSE)); sets
}
compute_Fcap_point <- function(model, df, sens_vars, times_for_prob, min_group_n = MIN_GROUP_N) {
  df2 <- complete_surv(df); if (nrow(df2) == 0L) return(NA_real_)
  keep <- stats::complete.cases(df2[, sens_vars, drop = FALSE])
  df2  <- df2[keep,,drop=FALSE]
  key  <- interaction(df2[, sens_vars, drop = FALSE], drop = TRUE, sep = "×")
  grp_index <- split(seq_len(nrow(df2)), key)
  S <- surv_mat(model, df2, times_for_prob)  # N x T
  Fcaps <- vapply(seq_along(times_for_prob), function(j){
    means <- vapply(grp_index, function(idx){
      if (length(idx) < min_group_n) return(NA_real_)
      mean(S[idx, j], na.rm = TRUE)
    }, numeric(1))
    means <- means[is.finite(means)]
    if (length(means) < 2) return(NA_real_)
    log(max(means) / min(means))
  }, numeric(1))
  mean(Fcaps, na.rm = TRUE)
}
intersectional_point_table <- function(df_test, models, combos, times_for_prob) {
  rows <- list()
  for (vars in combos) {
    tag <- paste(vars, collapse = "×")
    for (mname in names(models)) {
      val <- compute_Fcap_point(models[[mname]], df_test, vars, times_for_prob)
      rows[[length(rows)+1]] <- data.frame(Combo = tag, Method = mname, F_cap = round(val, 4), check.names = FALSE)
    }
  }
  dplyr::bind_rows(rows)
}
bootstrap_intersectional <- function(df_test, models, combos, times_for_prob, B = FAIR_BOOT_N, ncores = NCORES) {
  pb <- txtProgressBar(min = 0, max = 100, style = 3); on.exit(close(pb), add = TRUE)
  step <- max(1L, floor(B/100L)); ids <- split(seq_len(B), ceiling(seq_len(B)/step))
  all_res <- list(); done <- 0L
  for (chunk in ids) {
    chunk_res <- mclapply(chunk, function(b){
      idx  <- sample.int(nrow(df_test), replace = TRUE)
      boot <- df_test[idx,,drop=FALSE]
      intersectional_point_table(boot, models, combos, times_for_prob)
    }, mc.cores = ncores)
    all_res <- c(all_res, chunk_res); done <- done + length(chunk)
    setTxtProgressBar(pb, floor(100*done/B))
  }
  big <- dplyr::bind_rows(Map(function(df, rep_id) cbind(df, .rep = rep_id), all_res, seq_along(all_res)))
  big %>%
    dplyr::group_by(Combo, Method) %>%
    dplyr::summarise(
      Fcap_mean  = mean(F_cap, na.rm = TRUE),
      Fcap_lower = quantile(F_cap, 0.025, na.rm = TRUE),
      Fcap_upper = quantile(F_cap, 0.975, na.rm = TRUE),
      Fcap_se    = sd(F_cap, na.rm = TRUE),
      Fcap_n     = sum(is.finite(F_cap)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(`Fcap [95%CI]` = sprintf("%.4f [%.4f, %.4f]", Fcap_mean, Fcap_lower, Fcap_upper))
}

ALL_COMBOS  <- build_intersect_sets(FAIR_VARS, kmax = INTERSECT_MAX_K)  # only 2‑way
COMBOS_2    <- ALL_COMBOS
inter2_point<- intersectional_point_table(test_data, models_list, COMBOS_2, times_prob)
inter2_boot <- bootstrap_intersectional(test_data, models_list, COMBOS_2, times_prob,
                                        B = FAIR_BOOT_N, ncores = NCORES)

# --------------------------- Split subgroup sheets
Fair_tauInf <- fairness_combined %>% dplyr::filter(Horizon=="tau=∞")
Fair_tau24  <- fairness_combined %>% dplyr::filter(Horizon=="tau=24")
Fair_tau36  <- fairness_combined %>% dplyr::filter(Horizon=="tau=36")
Fair_tau60  <- fairness_combined %>% dplyr::filter(Horizon=="tau=60")

# --------------------------- WRITE EXCEL #2 (Fairness)
write_xlsx(
  list(
    Subgroup_tauInf           = Fair_tauInf,
    Subgroup_tau24            = Fair_tau24,
    Subgroup_tau36            = Fair_tau36,
    Subgroup_tau60            = Fair_tau60,
    Fairness_CI_FCI_FCG_Point = fairness_point_tbl %>%
      dplyr::transmute(Sensitive, Method, `CI(%)` = CI_percent, `F_CI` = F_CI, `F_CG` = F_CG),
    Fairness_CI_FCI_FCG_Boot  = fairness_boot_tbl %>%
      dplyr::select(Sensitive, Method,
                    `CI(%) [95%CI]`, CI_mean, CI_lower, CI_upper, CI_se, CI_n,
                    `F_CI [95%CI]`,  FCI_mean, FCI_lower, FCI_upper, FCI_se, FCI_n,
                    `F_CG [95%CI]`,  FCG_mean, FCG_lower, FCG_upper, FCG_se, FCG_n),
    Intersectional2_Point     = inter2_point,
    Intersectional2_Boot      = inter2_boot %>%
      dplyr::select(Combo, Method, `Fcap [95%CI]`, Fcap_mean, Fcap_lower, Fcap_upper, Fcap_se, Fcap_n)
  ),
  path = OUT_XLS_FAIR
)
log_i("✅ Saved Fairness tables → %s", OUT_XLS_FAIR)

# --------------------------- Save models & test (compressed)
saveRDS(list(cox_model = cox_model, rsf_model = rsf_model, bb_model = bb_model,
             bb_baseline = attr(bb_model,"bb_baseline"), tau_test_inf = tau_test_inf),
        file = MODELS_RDS, compress = "xz")
saveRDS(test_data, TEST_RDS, compress = "xz")
log_i("💾 Saved models → %s", MODELS_RDS)
log_i("💾 Saved test set → %s", TEST_RDS)

# --------------------------- Final housekeeping
gc(); cat("\nSession info:\n"); print(sessionInfo())
# ============================================================












# ============================================================
#Fairness @ Horizons — CI(%), F_CI(t), F_CG(t), F_cap(t) + 95% CIs
#          for τ = ∞, 24, 36, 60 and for 3 models.
#          Appends new sheets to OUT_XLS_FAIR.
# ============================================================

suppressPackageStartupMessages({ library(dplyr); library(writexl); library(parallel) })

# ---------- Helper: risk at a single time t (event probability by time t)
risk_at_time <- function(model, df, t){
  # RSF: use our safe survival adapter
  if (inherits(model, "rfsrc")) {
    S <- rsf_predict_surv_safe(model, df, times = t)
    return(as.numeric(1 - S[,1,drop=TRUE]))
  }
  # Cox / BlackBoost (mboost): try riskRegression first
  v <- try({
    pr <- riskRegression::predictRisk(model, newdata = df, times = t)
    as.numeric(pr[,1,drop=TRUE])
  }, silent = TRUE)
  if (!inherits(v, "try-error") && length(v) == nrow(df) && all(is.finite(v))) return(v)
  # Fallback via pec::predictSurvProb
  S2 <- try(pec::predictSurvProb(model, newdata = df, times = t), silent = TRUE)
  if (!inherits(S2, "try-error")) return(as.numeric(1 - S2[,1,drop=TRUE]))
  # Final fallback: scalar risk
  scalar_risk(model, df, tau = t)
}

# ---------- Helper: group-wise C-index and CI(%)
group_CF <- function(risk, time, status, grp){
  levs <- levels(factor(grp))
  sapply(levs, function(L){
    idx <- which(grp == L)
    if (length(idx) < 30 || length(unique(status[idx])) < 2) return(NA_real_)
    as.numeric(Hmisc::rcorr.cens(-risk[idx], Surv(time[idx], status[idx]))["C Index"])
  })
}
compute_CI_at_time <- function(model, df, sens_var, t){
  df2 <- complete_surv(df)
  if (nrow(df2) < 30) return(NA_real_)
  r <- risk_at_time(model, df2, t)
  ok <- is.finite(r)
  if (!any(ok)) return(NA_real_)
  CF <- group_CF(r[ok], df2$Survival.months[ok], df2$event_status[ok], factor(df2[[sens_var]][ok]))
  100 * max(abs(outer(CF, CF, `-`)), na.rm = TRUE)
}

# ---------- Bootstrap engine for (CI, FCI, FCG) at a single t
bootstrap_fairness_by_var_at_t <- function(models, df, fair_vars, t, B = FAIR_BOOT_N, ncores = NCORES){
  pb <- txtProgressBar(min = 0, max = 100, style = 3); on.exit(close(pb), add = TRUE)
  step <- max(1L, floor(B/100L)); ids <- split(seq_len(B), ceiling(seq_len(B)/step))
  res_all <- list(); done <- 0L
  for (chunk in ids){
    chunk_res <- mclapply(chunk, function(b){
      idx  <- sample.int(nrow(df), replace = TRUE)
      boot <- df[idx,,drop=FALSE]
      rows <- list()
      for (v in fair_vars) if (v %in% names(boot)) for (mname in names(models)) {
        ci_val <- compute_CI_at_time(models[[mname]], boot, v, t)
        ff     <- compute_FCI_FCG(models[[mname]],   boot, v, times = c(t), gamma = GAMMA_F)
        rows[[length(rows)+1]] <- data.frame(
          Sensitive = v, Method = mname,
          CI_percent = ci_val, F_CI = ff["FCI"], F_CG = ff["FCG"],
          stringsAsFactors = FALSE, check.names = FALSE
        )
      }
      dplyr::bind_rows(rows)
    }, mc.cores = ncores)
    res_all <- c(res_all, chunk_res); done <- done + length(chunk)
    setTxtProgressBar(pb, floor(100*done/B))
  }
  big <- dplyr::bind_rows(res_all)
  
  big %>%
    dplyr::group_by(Sensitive, Method) %>%
    dplyr::summarise(
      CI_mean   = mean(CI_percent, na.rm = TRUE),
      CI_lower  = quantile(CI_percent, 0.025, na.rm = TRUE),
      CI_upper  = quantile(CI_percent, 0.975, na.rm = TRUE),
      CI_se     = sd(CI_percent, na.rm = TRUE),
      CI_n      = sum(is.finite(CI_percent)),
      FCI_mean  = mean(F_CI, na.rm = TRUE),
      FCI_lower = quantile(F_CI, 0.025, na.rm = TRUE),
      FCI_upper = quantile(F_CI, 0.975, na.rm = TRUE),
      FCI_se    = sd(F_CI, na.rm = TRUE),
      FCI_n     = sum(is.finite(F_CI)),
      FCG_mean  = mean(F_CG, na.rm = TRUE),
      FCG_lower = quantile(F_CG, 0.025, na.rm = TRUE),
      FCG_upper = quantile(F_CG, 0.975, na.rm = TRUE),
      FCG_se    = sd(F_CG, na.rm = TRUE),
      FCG_n     = sum(is.finite(F_CG)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      `CI(%) [95%CI]` = sprintf("%.3f [%.3f, %.3f]", CI_mean,  CI_lower,  CI_upper),
      `F_CI [95%CI]`  = sprintf("%.4f [%.4f, %.4f]", FCI_mean, FCI_lower, FCI_upper),
      `F_CG [95%CI]`  = sprintf("%.4f [%.4f, %.4f]", FCG_mean, FCG_lower, FCG_upper)
    )
}

# ---------- Bootstrap for INTERSECTIONAL Fcap at a single t
bootstrap_intersectional_at_t <- function(models, df, combos, t, B = FAIR_BOOT_N, ncores = NCORES){
  pb <- txtProgressBar(min = 0, max = 100, style = 3); on.exit(close(pb), add = TRUE)
  step <- max(1L, floor(B/100L)); ids <- split(seq_len(B), ceiling(seq_len(B)/step))
  all_res <- list(); done <- 0L
  for (chunk in ids){
    chunk_res <- mclapply(chunk, function(b){
      idx  <- sample.int(nrow(df), replace = TRUE)
      boot <- df[idx,,drop=FALSE]
      intersectional_point_table(boot, models, combos, times_for_prob = c(t))
    }, mc.cores = ncores)
    all_res <- c(all_res, chunk_res); done <- done + length(chunk)
    setTxtProgressBar(pb, floor(100*done/B))
  }
  big <- dplyr::bind_rows(all_res)
  big %>%
    dplyr::group_by(Combo, Method) %>%
    dplyr::summarise(
      Fcap_mean  = mean(F_cap, na.rm = TRUE),
      Fcap_lower = quantile(F_cap, 0.025, na.rm = TRUE),
      Fcap_upper = quantile(F_cap, 0.975, na.rm = TRUE),
      Fcap_se    = sd(F_cap, na.rm = TRUE),
      Fcap_n     = sum(is.finite(F_cap)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(`Fcap [95%CI]` = sprintf("%.4f [%.4f, %.4f]", Fcap_mean, Fcap_lower, Fcap_upper))
}

# ---------- Orchestration across horizons
HZS_ALL <- c(Inf, HORIZONS)  # Inf, 24, 36, 60
tau_inf <- choose_tau_inf(test_data)
hz_to_time <- function(hz){
  if (is.infinite(hz)) return(tau_inf)
  min(choose_tau_inf(test_data), hz)
}

ALL_COMBOS_2 <- build_intersect_sets(FAIR_VARS, kmax = 2)

# Collect sheets to append
new_sheets <- list()

for (hz in HZS_ALL){
  t_use <- hz_to_time(hz)
  tag   <- if (is.infinite(hz)) "tauInf" else paste0("tau", hz)
  
  message(sprintf("[Fairness@%s] Bootstrap CI/FCI/FCG ...", tag))
  fair_boot_hz <- bootstrap_fairness_by_var_at_t(
    models = models_list, df = test_data, fair_vars = FAIR_VARS,
    t = t_use, B = FAIR_BOOT_N, ncores = NCORES
  ) %>%
    dplyr::select(
      Sensitive, Method,
      `CI(%) [95%CI]`, CI_mean, CI_lower, CI_upper, CI_se, CI_n,
      `F_CI [95%CI]`,  FCI_mean, FCI_lower, FCI_upper, FCI_se, FCI_n,
      `F_CG [95%CI]`,  FCG_mean, FCG_lower, FCG_upper, FCG_se, FCG_n
    )
  
  message(sprintf("[Fairness@%s] Bootstrap Intersectional Fcap ...", tag))
  inter_boot_hz <- bootstrap_intersectional_at_t(
    models = models_list, df = test_data, combos = ALL_COMBOS_2,
    t = t_use, B = FAIR_BOOT_N, ncores = NCORES
  ) %>%
    dplyr::select(Combo, Method, `Fcap [95%CI]`, Fcap_mean, Fcap_lower, Fcap_upper, Fcap_se, Fcap_n)
  
  # Add to sheet list with clear names
  new_sheets[[paste0("Fairness_Hz_", tag, "_CI_FCI_FCG")]] <- fair_boot_hz
  new_sheets[[paste0("Intersectional_Hz_", tag)]]           <- inter_boot_hz
}

# ---------- Append sheets to the existing OUT_XLS_FAIR
# (Read current workbook safely; if it doesn't exist, create new.)
append_to_excel <- function(path, new_tabs){
  old_tabs <- list()
  if (file.exists(path)){
    # Try to read; if fails, overwrite with new
    ok <- TRUE
    old_tabs <- tryCatch({
      readxl::excel_sheets(path) %>%
        setNames(.) %>%
        lapply(function(sh) readxl::read_xlsx(path, sheet = sh))
    }, error = function(e){ message("Could not read existing XLSX; will overwrite."); ok <<- FALSE; list() })
    if (!ok) old_tabs <- list()
  }
  out <- c(old_tabs, new_tabs)
  writexl::write_xlsx(out, path = path)
}

append_to_excel(OUT_XLS_FAIR, new_sheets)
log_i("📄 Appended horizon-specific fairness sheets → %s", OUT_XLS_FAIR)
# ============================================================

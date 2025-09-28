# =============================== Cox (multivariate) + 10-fold OOF + Median split + KM/log-rank ===============================
# Description:
# - Multivariate Cox model via glmnet (lasso) with 10-fold outer CV
# - Per-fold: fit on training fold (lambda chosen by cv.glmnet inside the fold), predict linear predictor on test fold
# - OOF linear predictor ("predicted hazard") -> median split (High vs Low risk)
# - Kaplan–Meier curves (log-rank test) based on OOF risk groups
# - Outputs: KM PDF/PNG, OOF risk CSV, log-rank result, per-fold selected features
#
# Put this script in your repo and run: Rscript survival_cv_pipeline.R
# Or source() it in an interactive R session.
# =============================================================================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(survival)
  library(survminer)
  library(glmnet)
  library(caret)
  library(ggplot2)
  library(readr)
  library(tibble)
})

set.seed(20250917)

# ------------------------------- Paths (relative; GitHub-friendly) -------------------------------
input_path <- file.path("data", "table3_merged.xlsx")  # <- put your file here
out_dir    <- file.path("outputs", "survival_cv")      # all results will be saved here
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- Columns -------------------------------
time_col  <- "follow_up_time_days"   # survival time (days)
event_col <- "END"                   # 0/1
id_col    <- NULL                    # optional: set to your subject ID column name if you have one (e.g., "ID")

# Feature set (edit to your variables)
features <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP","L.GC.ML.DG",
  "L.CA3","L.CA4","L.Whole_hippocampus","L.Hippocampal_tail","L.parasubiculum",
  "L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP","R.GC.ML.DG",
  "R.CA3","R.CA4","R.Whole_hippocampus","R.fimbria","R.parasubiculum",
  "R.Hippocampal_tail","R.HATA"
)

# ------------------------------- Read & clean -------------------------------
stopifnot(file.exists(input_path))
df <- as.data.frame(readxl::read_excel(input_path))
names(df) <- trimws(names(df))

# pick/create ID
if (!is.null(id_col) && id_col %in% names(df)) {
  df$.ID <- as.character(df[[id_col]])
} else if ("ID" %in% names(df)) {
  df$.ID <- as.character(df[["ID"]])
} else {
  df$.ID <- sprintf("row_%d", seq_len(nrow(df)))
}

# ensure time numeric and positive
df[[time_col]]  <- suppressWarnings(readr::parse_number(as.character(df[[time_col]])))
df[[time_col]][df[[time_col]] <= 0] <- NA_real_

# ensure event is strict 0/1 numeric
df[[event_col]] <- suppressWarnings(as.numeric(df[[event_col]]))
df[[event_col]][is.na(df[[event_col]])] <- 0
df[[event_col]][df[[event_col]] != 0]   <- 1

# keep only columns we need
keep_cols <- unique(c(".ID", time_col, event_col, features))
df <- df[, intersect(keep_cols, names(df)), drop = FALSE]

# simple median imputation for features (to avoid dropping many subjects)
is_feat <- intersect(features, names(df))
for (vn in is_feat) {
  if (!is.numeric(df[[vn]])) df[[vn]] <- suppressWarnings(as.numeric(df[[vn]]))
  if (anyNA(df[[vn]])) {
    med <- suppressWarnings(median(df[[vn]], na.rm = TRUE))
    if (is.finite(med)) df[[vn]][is.na(df[[vn]])] <- med
  }
}

# require non-missing time/event and at least some features
df <- df[!is.na(df[[time_col]]) & !is.na(df[[event_col]]), , drop = FALSE]

# if still missing features after imputation, drop incomplete rows for modeling matrix
df <- df[stats::complete.cases(df[, is_feat, drop = FALSE]), , drop = FALSE]

if (nrow(df) < 30 || sum(df[[event_col]] == 1) < 5) {
  stop("Not enough samples or events for CV Cox modeling (need roughly >=30 rows and >=5 events).")
}

# ------------------------------- 10-fold CV to get OOF linear predictor -------------------------------
folds <- caret::createFolds(df$.ID, k = 10, returnTrain = TRUE)
oof_lp <- rep(NA_real_, nrow(df))
sel_feats_per_fold <- vector("list", length(folds))

for (i in seq_along(folds)) {
  tr_idx <- folds[[i]]
  te_idx <- setdiff(seq_len(nrow(df)), tr_idx)
  dtr <- df[tr_idx, , drop = FALSE]
  dte <- df[te_idx, , drop = FALSE]
  
  # glmnet cox design
  x_tr <- model.matrix(~ . - 1,
                       data = dtr[, is_feat, drop = FALSE])  # -1 no intercept
  y_tr <- with(dtr, Surv(get(time_col), get(event_col)))
  
  # inner CV on the training fold to pick lambda (no leakage)
  cvfit <- glmnet::cv.glmnet(
    x = x_tr, y = y_tr, family = "cox",
    nfolds = 5, alpha = 1, standardize = TRUE
  )
  
  # training a Cox model implicitly done by cv.glmnet; we just use coef at lambda.min
  beta <- as.matrix(coef(cvfit, s = "lambda.min"))
  keep <- rownames(beta)[as.numeric(beta) != 0]
  sel_feats_per_fold[[i]] <- setdiff(keep, "(Intercept)")
  
  # predict LP on test fold (OOF)
  x_te <- model.matrix(~ . - 1,
                       data = dte[, is_feat, drop = FALSE])
  lp_te <- as.numeric(predict(cvfit, newx = x_te, s = "lambda.min", type = "link"))
  oof_lp[te_idx] <- lp_te
}

# sanity check
if (sum(is.finite(oof_lp)) < 0.8 * nrow(df)) {
  warning("A lot of OOF predictions are NA; check folds or data quality.")
}

# ------------------------------- Dichotomize by OOF median & KM/log-rank -------------------------------
# high risk = lp > median(oof_lp)
med_lp <- stats::median(oof_lp, na.rm = TRUE)
risk_grp <- ifelse(oof_lp > med_lp, "High risk", "Low risk")
risk_grp[!is.finite(oof_lp)] <- NA_character_

dat_km <- df %>%
  transmute(
    ID = .ID,
    time = .data[[time_col]],
    event = .data[[event_col]],
    OOF_LP = oof_lp,
    RiskGroup = factor(risk_grp, levels = c("Low risk","High risk"))
  ) %>%
  filter(!is.na(RiskGroup))

# need both groups and some events
if (length(unique(dat_km$RiskGroup)) < 2 || sum(dat_km$event == 1) < 5) {
  stop("Not enough separation/events to draw KM with median split of OOF risk.")
}

sf <- survfit(Surv(time, event) ~ RiskGroup, data = dat_km)
lr <- survdiff(Surv(time, event) ~ RiskGroup, data = dat_km)
p_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

# ------------------------------- Save outputs -------------------------------
# 1) OOF risk CSV
readr::write_csv(dat_km, file.path(out_dir, "oof_risk_groups.csv"))

# 2) per-fold selected features
sink(file.path(out_dir, "selected_features_by_fold.txt"))
cat("Selected (non-zero) features by fold (lambda.min):\n\n")
for (i in seq_along(sel_feats_per_fold)) {
  cat(sprintf("Fold %02d: %s\n", i,
              if (length(sel_feats_per_fold[[i]]) == 0) "(none)"
              else paste(sel_feats_per_fold[[i]], collapse = ", ")))
}
sink()

# 3) KM plot (PDF + PNG)
p_km <- ggsurvplot(
  sf, data = dat_km,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, risk.table.height = 0.22,
  surv.median.line = "hv",
  legend.title = NULL, legend.labs = c("Low risk","High risk"),
  palette = c("#2C7BB6","#D7191C")
)
ggsave(file.path(out_dir, "KM_OOF_median_split.pdf"), p_km$plot, width = 7.5, height = 6.0, device = "pdf")
ggsave(file.path(out_dir, "KM_OOF_median_split.png"), p_km$plot, width = 7.5, height = 6.0, dpi = 300)

# 4) Log-rank summary
write_lines(c(
  "Cox (glmnet) multivariate with 10-fold outer CV",
  "OOF linear predictor used for median split (High vs Low risk)",
  sprintf("Log-rank p = %.4g", p_lr),
  sprintf("N = %d (events = %d)", nrow(dat_km), sum(dat_km$event == 1))
), file.path(out_dir, "logrank_summary.txt"))

# 5) Optional: OOF concordance (Harrell's C) on OOF LP
# Note: survConcordance expects a numeric predictor; larger = higher risk
c_idx <- tryCatch({
  sc <- survival::survConcordance(Surv(time, event) ~ OOF_LP, data = dat_km)
  sc$concordance
}, error = function(e) NA_real_)
write_lines(sprintf("OOF Harrell's C (approx): %.3f", c_idx),
            file.path(out_dir, "oof_concordance.txt"))

cat("\nDone. Results saved to: ", normalizePath(out_dir), "\n", sep = "")

# ============================================================
# 08_Progression_analysis.R
# Leakage-safe prognostic analysis for HC/MCI and MS endpoints:
# training-fold LASSO Cox, out-of-fold risk, C-index, median
# risk stratification, Kaplan-Meier curves, and log-rank tests.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Usage:
# Rscript 08_Progression_analysis.R <HC_MCI_data.csv|xlsx|-> \
#   <MS_data.csv|xlsx|-> <output_dir> [HC_MCI_features|auto] \
#   [MS_features] [seed]
#
# Defaults inherited from the previous scripts:
#   HC/MCI: time = follow_up_time_days; event = END;
#           features = auto-detected 24 hippocampal centiles.
#   MS:     time = Followup_time_month;
#           events = EDSS_Progress and SPMS_conversion;
#           features = EX.
#
# Comma-separate feature names. Use "-" to skip either dataset.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 08_Progression_analysis.R <HC_MCI_data.csv|xlsx|-> ",
    "<MS_data.csv|xlsx|-> <output_dir> [HC_MCI_features|auto] ",
    "[MS_features] [seed]"
  )
}

hc_mci_arg <- args[[1]]
ms_arg <- args[[2]]
output_dir <- args[[3]]
hc_mci_feature_arg <- if (length(args) >= 4) args[[4]] else "auto"
ms_feature_arg <- if (length(args) >= 5) args[[5]] else "EX"
random_seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20250917L

if (identical(hc_mci_arg, "-") && identical(ms_arg, "-")) {
  stop("At least one prognostic dataset must be supplied")
}
if (!is.finite(random_seed)) stop("seed must be an integer")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

required_packages <- c("survival", "glmnet", "caret")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

outer_folds <- 10L
inner_folds <- 5L
minimum_n <- 30L
minimum_events <- 5L
significance_threshold <- 0.05
set.seed(random_seed)

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

sanitize_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

normalize_name <- function(x) gsub("[^a-z0-9]", "", tolower(safe_chr(x)))

bind_rows_base <- function(x) {
  x <- Filter(function(item) !is.null(item) && nrow(item) > 0, x)
  if (length(x) == 0) data.frame() else do.call(rbind, x)
}

read_table_auto <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  extension <- tolower(tools::file_ext(path))
  if (extension == "csv") {
    data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (extension %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required for Excel prognostic data")
    }
    data <- as.data.frame(readxl::read_excel(path, sheet = 1), check.names = FALSE)
  } else {
    stop("Prognostic data must be CSV, XLSX, or XLS: ", path)
  }
  names(data) <- trimws(names(data))
  attr(data, "source_path") <- path
  data
}

resolve_column <- function(data, requested) {
  if (requested %in% names(data)) return(requested)
  match_index <- match(normalize_name(requested), normalize_name(names(data)))
  if (is.na(match_index)) stop("Required column not found: ", requested)
  names(data)[match_index]
}

old_hippocampal_features <- c(
  "L.subiculum", "L.CA1", "L.presubiculum", "L.molecular_layer_HP",
  "L.GC.ML.DG", "L.CA3", "L.CA4", "L.Whole_hippocampus",
  "L.Hippocampal_tail", "L.parasubiculum", "L.fimbria", "L.HATA",
  "R.subiculum", "R.CA1", "R.presubiculum", "R.molecular_layer_HP",
  "R.GC.ML.DG", "R.CA3", "R.CA4", "R.Whole_hippocampus",
  "R.fimbria", "R.parasubiculum", "R.Hippocampal_tail", "R.HATA"
)

new_hippocampal_features <- unlist(
  lapply(c("left", "right"), function(side) paste0(side, "_", c(
    "Whole_hippocampus", "subiculum", "CA1", "presubiculum",
    "molecular_layer_HP", "GC.ML.DG", "CA3", "CA4",
    "Hippocampal_tail", "parasubiculum", "fimbria", "HATA"
  ))),
  use.names = FALSE
)

resolve_features <- function(data, feature_arg, task_type) {
  if (tolower(feature_arg) != "auto") {
    requested <- trimws(unlist(strsplit(feature_arg, ",", fixed = TRUE)))
    requested <- requested[requested != ""]
    resolved <- vapply(requested, function(feature) {
      if (feature %in% names(data)) return(feature)
      index <- match(normalize_name(feature), normalize_name(names(data)))
      if (is.na(index)) stop("Feature not found: ", feature)
      names(data)[index]
    }, character(1))
    return(unique(resolved))
  }

  old_present <- old_hippocampal_features[old_hippocampal_features %in% names(data)]
  if (length(old_present) == length(old_hippocampal_features)) return(old_present)
  new_present <- new_hippocampal_features[new_hippocampal_features %in% names(data)]
  if (length(new_present) == length(new_hippocampal_features)) return(new_present)
  centile_columns <- grep("_centile$", names(data), value = TRUE, ignore.case = TRUE)
  if (length(centile_columns) > 0) return(centile_columns)
  stop("Could not auto-detect prognostic features for task: ", task_type)
}

choose_id <- function(data) {
  candidates <- c("match_id", "ID", ".ID", "Subject_ID", "subject_id")
  existing <- candidates[candidates %in% names(data)]
  if (length(existing) > 0) {
    id <- safe_chr(data[[existing[[1]]]])
    id[id == ""] <- paste0("row_", which(id == ""))
    return(id)
  }
  paste0("row_", seq_len(nrow(data)))
}

make_outer_test_folds <- function(id, event) {
  subject_event <- tapply(event, id, max, na.rm = TRUE)
  set.seed(random_seed)
  subject_folds <- caret::createFolds(
    factor(subject_event, levels = c(0, 1)),
    k = outer_folds,
    returnTrain = FALSE
  )
  lapply(subject_folds, function(subject_index) {
    which(id %in% names(subject_event)[subject_index])
  })
}

prepare_fold_matrices <- function(training_data, test_data, features) {
  x_train <- as.data.frame(
    training_data[, features, drop = FALSE], check.names = FALSE
  )
  x_test <- as.data.frame(
    test_data[, features, drop = FALSE], check.names = FALSE
  )
  if (any(!is.finite(as.matrix(x_train))) || any(!is.finite(as.matrix(x_test)))) {
    stop("Non-finite predictors remained after complete-case filtering")
  }

  nonconstant <- features[vapply(x_train, function(x) {
    is.finite(stats::sd(x)) && stats::sd(x) > 0
  }, logical(1))]
  if (length(nonconstant) == 0) return(NULL)
  list(
    train = as.matrix(x_train[, nonconstant, drop = FALSE]),
    test = as.matrix(x_test[, nonconstant, drop = FALSE]),
    features = nonconstant
  )
}

draw_km_plot <- function(fit, logrank_p, task_label, time_label, output_path, device) {
  if (device == "png") {
    grDevices::png(output_path, width = 2100, height = 1650, res = 300)
  } else {
    grDevices::pdf(output_path, width = 7, height = 5.5)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    fit,
    col = c("#2C7BB6", "#D7191C"),
    lwd = 2,
    xlab = time_label,
    ylab = "Survival probability",
    mark.time = TRUE,
    conf.int = FALSE,
    main = task_label
  )
  graphics::legend(
    "bottomleft",
    legend = c("Low risk", "High risk"),
    col = c("#2C7BB6", "#D7191C"),
    lwd = 2,
    bty = "n"
  )
  graphics::mtext(
    sprintf("Two-sided log-rank P = %.4g", logrank_p),
    side = 3,
    line = 0.2,
    cex = 0.85
  )
}

run_endpoint <- function(data, dataset_name, time_requested, event_requested, features) {
  source_path <- attr(data, "source_path")
  time_column <- resolve_column(data, time_requested)
  event_column <- resolve_column(data, event_requested)
  task_name <- paste(dataset_name, event_requested, sep = "__")
  task_stub <- sanitize_name(task_name)
  message("Prognostic task: ", task_name)

  d <- data.frame(
    ID = choose_id(data),
    time = safe_num(data[[time_column]]),
    event = safe_num(data[[event_column]]),
    data[, features, drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (feature in features) d[[feature]] <- safe_num(d[[feature]])
  d$event[is.finite(d$event) & d$event != 0] <- 1

  # Apply complete-case eligibility before constructing any CV folds:
  # first require a valid survival outcome, then require every model
  # predictor to be finite. No predictor imputation is performed.
  d <- d[
    is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1),
    , drop = FALSE
  ]
  complete_predictors <- rowSums(is.finite(as.matrix(
    d[, features, drop = FALSE]
  ))) == length(features)
  d <- d[complete_predictors, , drop = FALSE]

  n_subjects <- length(unique(d$ID))
  n_events <- sum(d$event == 1)
  if (nrow(d) < minimum_n || n_subjects < outer_folds || n_events < minimum_events) {
    return(list(
      status = data.frame(
        task = task_name, status = "skipped",
        reason = "insufficient observations, subjects, or events",
        n = nrow(d), n_subjects = n_subjects, events = n_events,
        stringsAsFactors = FALSE
      )
    ))
  }

  test_folds <- make_outer_test_folds(d$ID, d$event)
  oof_risk <- rep(NA_real_, nrow(d))
  fold_rows <- list()
  coefficient_rows <- list()

  for (fold_index in seq_along(test_folds)) {
    test_index <- test_folds[[fold_index]]
    training_index <- setdiff(seq_len(nrow(d)), test_index)
    training_data <- d[training_index, , drop = FALSE]
    test_data <- d[test_index, , drop = FALSE]
    training_events <- sum(training_data$event == 1)
    test_events <- sum(test_data$event == 1)

    fold_status <- "completed"
    fold_reason <- ""
    selected_features <- character(0)
    lambda_min <- NA_real_
    matrices <- prepare_fold_matrices(training_data, test_data, features)
    if (is.null(matrices) || training_events < 2) {
      fold_status <- "failed"
      fold_reason <- "no usable training features or fewer than two training events"
    } else {
      survival_training <- survival::Surv(training_data$time, training_data$event)
      set.seed(random_seed + fold_index)
      inner_fold_id <- as.integer(caret::createFolds(
        factor(training_data$event, levels = c(0, 1)),
        k = inner_folds,
        list = FALSE
      ))
      fit_error <- NULL
      fit <- tryCatch(
        glmnet::cv.glmnet(
          x = matrices$train,
          y = survival_training,
          family = "cox",
          alpha = 1,
          foldid = inner_fold_id,
          standardize = TRUE,
          type.measure = "deviance"
        ),
        error = function(error) {
          fit_error <<- conditionMessage(error)
          NULL
        }
      )
      if (is.null(fit)) {
        fold_status <- "failed"
        fold_reason <- paste0("training-fold LASSO Cox failed: ", fit_error)
      } else {
        lambda_min <- fit$lambda.min
        coefficients <- as.matrix(stats::coef(fit, s = "lambda.min"))[, 1]
        selected_features <- names(coefficients)[
          is.finite(coefficients) & abs(coefficients) > 0
        ]
        if (length(selected_features) == 0) {
          fold_status <- "failed"
          fold_reason <- "LASSO selected no features at lambda.min"
        } else {
          test_risk <- tryCatch(
            as.numeric(stats::predict(
              fit,
              newx = matrices$test,
              s = "lambda.min",
              type = "link"
            )),
            error = function(error) NULL
          )
          if (is.null(test_risk) || length(test_risk) != length(test_index) ||
              any(!is.finite(test_risk))) {
            fold_status <- "failed"
            fold_reason <- "training-fold Cox coefficients could not score the held-out fold"
          } else {
            oof_risk[test_index] <- test_risk
            coefficient_rows[[length(coefficient_rows) + 1]] <- data.frame(
              task = task_name,
              outer_fold = fold_index,
              feature = selected_features,
              coefficient = as.numeric(coefficients[selected_features]),
              lambda_min = lambda_min,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    fold_rows[[length(fold_rows) + 1]] <- data.frame(
      task = task_name,
      outer_fold = fold_index,
      training_n = length(training_index),
      training_events = training_events,
      heldout_n = length(test_index),
      heldout_events = test_events,
      n_usable_training_features = if (is.null(matrices)) 0L else length(matrices$features),
      n_selected_features = length(selected_features),
      selected_features = paste(selected_features, collapse = ";"),
      lambda_min = lambda_min,
      status = fold_status,
      reason = fold_reason,
      stringsAsFactors = FALSE
    )
  }

  fold_audit <- bind_rows_base(fold_rows)
  coefficients <- bind_rows_base(coefficient_rows)
  if (any(!is.finite(oof_risk))) {
    return(list(
      status = data.frame(
        task = task_name,
        status = "failed",
        reason = "one or more outer folds failed; no full-sample model was used to fill OOF risk",
        n = nrow(d), n_subjects = n_subjects, events = n_events,
        stringsAsFactors = FALSE
      ),
      fold_audit = fold_audit,
      coefficients = coefficients
    ))
  }

  median_risk <- stats::median(oof_risk)
  risk_group <- factor(
    ifelse(oof_risk > median_risk, "High risk", "Low risk"),
    levels = c("Low risk", "High risk")
  )
  if (nlevels(droplevels(risk_group)) < 2) {
    return(list(
      status = data.frame(
        task = task_name, status = "failed",
        reason = "median OOF risk split produced fewer than two groups",
        n = nrow(d), n_subjects = n_subjects, events = n_events,
        stringsAsFactors = FALSE
      ),
      fold_audit = fold_audit,
      coefficients = coefficients
    ))
  }

  oof_data <- data.frame(
    task = task_name,
    ID = d$ID,
    time = d$time,
    event = d$event,
    OOF_risk_score = oof_risk,
    OOF_median_risk = median_risk,
    Risk_group = risk_group,
    stringsAsFactors = FALSE
  )
  concordance_fit <- survival::concordance(
    survival::Surv(time, event) ~ OOF_risk_score,
    data = oof_data,
    reverse = TRUE
  )
  c_index <- unname(concordance_fit$concordance)
  c_se <- if (!is.null(concordance_fit$std.err)) {
    unname(concordance_fit$std.err)
  } else if (!is.null(concordance_fit$var)) {
    sqrt(unname(concordance_fit$var))
  } else {
    NA_real_
  }

  km_fit <- survival::survfit(
    survival::Surv(time, event) ~ Risk_group,
    data = oof_data
  )
  logrank <- survival::survdiff(
    survival::Surv(time, event) ~ Risk_group,
    data = oof_data
  )
  logrank_df <- length(logrank$n) - 1
  logrank_p <- stats::pchisq(logrank$chisq, df = logrank_df, lower.tail = FALSE)
  km_summary <- summary(km_fit)
  km_table <- data.frame(
    task = task_name,
    stratum = if (is.null(km_summary$strata)) "All" else as.character(km_summary$strata),
    time = km_summary$time,
    n_risk = km_summary$n.risk,
    n_event = km_summary$n.event,
    survival = km_summary$surv,
    standard_error = km_summary$std.err,
    lower_95CI = km_summary$lower,
    upper_95CI = km_summary$upper,
    stringsAsFactors = FALSE
  )

  task_dir <- file.path(output_dir, task_stub)
  if (!dir.exists(task_dir)) dir.create(task_dir, recursive = TRUE)
  write.csv(oof_data, file.path(task_dir, "OOF_risk_scores.csv"), row.names = FALSE)
  write.csv(fold_audit, file.path(task_dir, "training_fold_audit.csv"), row.names = FALSE)
  write.csv(coefficients, file.path(task_dir, "selected_features_and_Cox_coefficients.csv"), row.names = FALSE)
  write.csv(km_table, file.path(task_dir, "Kaplan_Meier_estimates.csv"), row.names = FALSE)
  draw_km_plot(
    km_fit, logrank_p, task_name, time_column,
    file.path(task_dir, "Kaplan_Meier_OOF_median_risk.png"), "png"
  )
  draw_km_plot(
    km_fit, logrank_p, task_name, time_column,
    file.path(task_dir, "Kaplan_Meier_OOF_median_risk.pdf"), "pdf"
  )

  performance <- data.frame(
    task = task_name,
    dataset = dataset_name,
    source_file = source_path,
    time_column = time_column,
    event_column = event_column,
    n = nrow(oof_data),
    n_subjects = n_subjects,
    events = sum(oof_data$event == 1),
    n_features_entered = length(features),
    OOF_C_index = c_index,
    OOF_C_index_95CI_low = if (is.finite(c_se)) max(0, c_index - 1.96 * c_se) else NA_real_,
    OOF_C_index_95CI_high = if (is.finite(c_se)) min(1, c_index + 1.96 * c_se) else NA_real_,
    OOF_median_risk = median_risk,
    logrank_chisq = unname(logrank$chisq),
    logrank_df = logrank_df,
    logrank_p_two_sided = logrank_p,
    significant_at_0_05 = is.finite(logrank_p) && logrank_p < significance_threshold,
    stringsAsFactors = FALSE
  )
  write.csv(performance, file.path(task_dir, "prognostic_performance.csv"), row.names = FALSE)
  list(
    status = data.frame(
      task = task_name, status = "completed", reason = "",
      n = nrow(d), n_subjects = n_subjects, events = n_events,
      stringsAsFactors = FALSE
    ),
    performance = performance,
    fold_audit = fold_audit,
    coefficients = coefficients
  )
}

all_status <- list()
all_performance <- list()
all_fold_audit <- list()
all_coefficients <- list()
settings_rows <- list()

append_result <- function(result) {
  if (!is.null(result$status)) all_status[[length(all_status) + 1]] <<- result$status
  if (!is.null(result$performance)) all_performance[[length(all_performance) + 1]] <<- result$performance
  if (!is.null(result$fold_audit)) all_fold_audit[[length(all_fold_audit) + 1]] <<- result$fold_audit
  if (!is.null(result$coefficients)) all_coefficients[[length(all_coefficients) + 1]] <<- result$coefficients
}

if (!identical(hc_mci_arg, "-")) {
  hc_mci <- read_table_auto(hc_mci_arg)
  hc_features <- resolve_features(hc_mci, hc_mci_feature_arg, "HC_MCI")
  append_result(run_endpoint(
    hc_mci,
    dataset_name = "HC_MCI",
    time_requested = "follow_up_time_days",
    event_requested = "END",
    features = hc_features
  ))
  settings_rows[[length(settings_rows) + 1]] <- data.frame(
    dataset = "HC_MCI", source_file = attr(hc_mci, "source_path"),
    time_column = "follow_up_time_days", event_columns = "END",
    features = paste(hc_features, collapse = ";"), stringsAsFactors = FALSE
  )
}

if (!identical(ms_arg, "-")) {
  ms <- read_table_auto(ms_arg)
  ms_features <- resolve_features(ms, ms_feature_arg, "MS")
  for (endpoint in c("EDSS_Progress", "SPMS_conversion")) {
    append_result(run_endpoint(
      ms,
      dataset_name = "MS",
      time_requested = "Followup_time_month",
      event_requested = endpoint,
      features = ms_features
    ))
  }
  settings_rows[[length(settings_rows) + 1]] <- data.frame(
    dataset = "MS", source_file = attr(ms, "source_path"),
    time_column = "Followup_time_month",
    event_columns = "EDSS_Progress;SPMS_conversion",
    features = paste(ms_features, collapse = ";"), stringsAsFactors = FALSE
  )
}

status_table <- bind_rows_base(all_status)
performance_table <- bind_rows_base(all_performance)
fold_audit_table <- bind_rows_base(all_fold_audit)
coefficient_table <- bind_rows_base(all_coefficients)
dataset_settings <- bind_rows_base(settings_rows)

if (nrow(status_table) > 0) {
  write.csv(status_table, file.path(output_dir, "progression_task_status.csv"), row.names = FALSE)
}
if (nrow(performance_table) > 0) {
  write.csv(
    performance_table,
    file.path(output_dir, "progression_OOF_prognostic_performance.csv"),
    row.names = FALSE
  )
}
if (nrow(fold_audit_table) > 0) {
  write.csv(
    fold_audit_table,
    file.path(output_dir, "progression_all_training_fold_audit.csv"),
    row.names = FALSE
  )
}
if (nrow(coefficient_table) > 0) {
  write.csv(
    coefficient_table,
    file.path(output_dir, "progression_all_selected_features_Cox_coefficients.csv"),
    row.names = FALSE
  )
}

method_settings <- data.frame(
  item = c(
    "outer cross-validation", "inner LASSO tuning", "missing-data handling",
    "feature standardization", "OOF risk", "failed outer folds",
    "prognostic performance", "risk stratification", "survival visualization",
    "survival comparison", "statistical tests", "significance threshold", "seed"
  ),
  value = c(
    "10-fold; repeated rows from one ID remain in the same held-out fold",
    "5-fold cv.glmnet within each outer training fold; alpha=1; lambda.min",
    paste0(
      "complete-case filtering before cross-validation: valid time/event and all ",
      "model predictors required; no imputation"
    ),
    "glmnet standardization estimated from each training fold only",
    "training-fold LASSO Cox coefficients applied only to the corresponding held-out fold",
    "recorded as failures; never filled using a full-sample model",
    "Harrell concordance index from pooled out-of-fold risk scores",
    "median pooled out-of-fold risk score",
    "Kaplan-Meier curves for low- and high-risk groups",
    "two-sided log-rank test",
    "two-sided unless otherwise specified",
    significance_threshold,
    random_seed
  ),
  stringsAsFactors = FALSE
)
write.csv(method_settings, file.path(output_dir, "progression_method_settings.csv"), row.names = FALSE)
if (nrow(dataset_settings) > 0) {
  write.csv(dataset_settings, file.path(output_dir, "progression_dataset_settings.csv"), row.names = FALSE)
}

message(
  "Progression analysis completed for ", nrow(performance_table),
  " endpoint(s). Outputs: ", output_dir
)

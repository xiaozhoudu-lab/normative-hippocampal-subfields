# ============================================================
# 07_Cognitive_prediction.R
# Part 1. Cognitive associations
# Part 2. Cross-sectional cognitive prediction
# Part 3. Longitudinal cognitive prediction
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Usage:
# Rscript 07_Cognitive_prediction.R <normative_scores_long.csv> \
#   <cognitive_data.csv|xlsx> <output_dir> [cognitive_outcomes] \
#   [follow_up_time_column] [diagnosis_column] [healthy_control_labels] \
#   [prediction_groups|auto] [seed] [individual_CMD.csv|-]
#
# Defaults:
#   cognitive_outcomes       = MOCA,MMSE,SDMT,BVMT,PASAT
#   follow_up_time_column    = follow_up_time_days
#   diagnosis_column         = Diagnosis
#   healthy_control_labels   = HC
#   prediction_groups        = auto (all non-HC diagnosis groups)
#   individual_CMD.csv       = - (regional associations only)
#
# Prediction method:
#   - radial-basis-function support vector regression;
#   - explicit outer 10-fold cross-validation;
#   - inner 5-fold tuning using each outer training set only;
#   - centering/scaling estimated inside training resamples only;
#   - age, sex, and 24 hippocampal centile features as predictors;
#   - pooled held-out predictions evaluated with Pearson correlation.
# Longitudinal folds are grouped by match_id to prevent subject leakage.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 07_Cognitive_prediction.R <normative_scores_long.csv> ",
    "<cognitive_data.csv|xlsx> <output_dir> [cognitive_outcomes] ",
    "[follow_up_time_column] [diagnosis_column] [healthy_control_labels] ",
    "[prediction_groups|auto] [seed] [individual_CMD.csv|-]"
  )
}

score_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
cognitive_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
outcome_arg <- if (length(args) >= 4) args[[4]] else "MOCA,MMSE,SDMT,BVMT,PASAT"
follow_up_time_column <- if (length(args) >= 5) args[[5]] else "follow_up_time_days"
diagnosis_column <- if (length(args) >= 6) args[[6]] else "Diagnosis"
healthy_arg <- if (length(args) >= 7) args[[7]] else "HC"
prediction_group_arg <- if (length(args) >= 8) args[[8]] else "auto"
random_seed <- if (length(args) >= 9) as.integer(args[[9]]) else 20260527L
cmd_arg <- if (length(args) >= 10) args[[10]] else "-"

if (!is.finite(random_seed)) stop("seed must be an integer")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

required_packages <- c("caret", "kernlab")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

outer_folds <- 10L
inner_folds <- 5L
minimum_association_n <- 10L
minimum_prediction_n <- 30L
set.seed(random_seed)

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

make_match_id <- function(df) {
  if ("match_id" %in% names(df)) return(safe_chr(df$match_id))
  candidates <- list(
    c("Freesufer_Path2", "Freesufer_Path3"),
    c("Freesurfer_Path2", "Freesurfer_Path3")
  )
  for (pair in candidates) {
    if (all(pair %in% names(df))) {
      return(paste0(safe_chr(df[[pair[[1]]]]), safe_chr(df[[pair[[2]]]])))
    }
  }
  stop(
    "Cognitive data require either match_id or the same FreeSurfer path columns ",
    "used by 02_Apply_normative_model.R"
  )
}

normalize_diagnosis <- function(x) {
  out <- safe_chr(x)
  out[toupper(out) == "SVD"] <- "CSVD"
  out
}

normalize_hemisphere <- function(x) {
  key <- gsub("[^a-z]", "", tolower(safe_chr(x)))
  out <- rep(NA_character_, length(key))
  out[key %in% c("l", "lh", "left")] <- "left"
  out[key %in% c("r", "rh", "right")] <- "right"
  out
}

normalize_roi <- function(x) {
  key <- gsub("[^a-z0-9]", "", tolower(safe_chr(x)))
  aliases <- c(
    wholehippocampus = "Whole_hippocampus",
    subiculum = "subiculum",
    ca1 = "CA1",
    presubiculum = "presubiculum",
    molecularlayerhp = "molecular_layer_HP",
    molecularlayer = "molecular_layer_HP",
    gcmldg = "GC.ML.DG",
    ca3 = "CA3",
    ca4 = "CA4",
    hippocampaltail = "Hippocampal_tail",
    parasubiculum = "parasubiculum",
    fimbria = "fimbria",
    hata = "HATA"
  )
  unname(aliases[key])
}

collapse_unique <- function(x) {
  values <- sort(unique(safe_chr(x)))
  values <- values[values != ""]
  if (length(values) == 0) "" else paste(values, collapse = ";")
}

bind_rows_base <- function(x) {
  x <- Filter(function(item) !is.null(item) && nrow(item) > 0, x)
  if (length(x) == 0) data.frame() else do.call(rbind, x)
}

read_table_auto <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension == "csv") {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (extension %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required for Excel cognitive data")
    }
    as.data.frame(readxl::read_excel(path, sheet = 1), check.names = FALSE)
  } else {
    stop("Cognitive data must be CSV, XLSX, or XLS")
  }
}

roi_order <- c(
  "Whole_hippocampus", "subiculum", "CA1", "presubiculum",
  "molecular_layer_HP", "GC.ML.DG", "CA3", "CA4",
  "Hippocampal_tail", "parasubiculum", "fimbria", "HATA"
)
feature_order <- unlist(
  lapply(c("left", "right"), function(side) paste0(side, "_", roi_order)),
  use.names = FALSE
)

# ------------------------------------------------------------
# Prepare 24 regional percentile features from script 02 output.
# ------------------------------------------------------------
scores <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
required_score_columns <- c("match_id", "hemisphere", "ROI", "Quant_score")
missing_score_columns <- setdiff(required_score_columns, names(scores))
if (length(missing_score_columns) > 0) {
  stop(
    "The score file must be normative_scores_long.csv from ",
    "02_Apply_normative_model.R; missing: ", paste(missing_score_columns, collapse = ", ")
  )
}
scores$match_id <- safe_chr(scores$match_id)
scores$hemisphere_cognition <- normalize_hemisphere(scores$hemisphere)
scores$ROI_cognition <- normalize_roi(scores$ROI)
if (any(scores$match_id == "")) stop("The score file contains empty match_id values")
if (anyNA(scores$hemisphere_cognition) || anyNA(scores$ROI_cognition)) {
  stop("The score file contains unrecognized hemisphere or ROI values")
}
scores$feature_name_cognition <- paste0(
  scores$hemisphere_cognition, "_", scores$ROI_cognition
)
unexpected_features <- setdiff(unique(scores$feature_name_cognition), feature_order)
if (length(unexpected_features) > 0) {
  stop("Unexpected cognitive features: ", paste(unexpected_features, collapse = ", "))
}
duplicate_key <- duplicated(scores[c("match_id", "feature_name_cognition")]) |
  duplicated(scores[c("match_id", "feature_name_cognition")], fromLast = TRUE)
if (any(duplicate_key)) stop("Duplicate participant-feature rows found in score input")

percentile <- safe_num(scores$Quant_score)
finite_percentile <- percentile[is.finite(percentile)]
if (length(finite_percentile) == 0) stop("Quant_score contains no finite values")
percentile_scale <- if (max(finite_percentile) > 1.5) "0_to_100" else "0_to_1"
if (percentile_scale == "0_to_100") percentile <- percentile / 100
if (any(is.finite(percentile) & (percentile < 0 | percentile > 1))) {
  stop("Quant_score contains values outside the detected percentile scale")
}

participant_ids <- unique(scores$match_id)
feature_matrix <- matrix(
  NA_real_, nrow = length(participant_ids), ncol = length(feature_order),
  dimnames = list(participant_ids, feature_order)
)
feature_matrix[cbind(
  match(scores$match_id, participant_ids),
  match(scores$feature_name_cognition, feature_order)
)] <- percentile

feature_data <- data.frame(
  match_id = participant_ids,
  n_hippocampal_features = rowSums(is.finite(feature_matrix)),
  complete_24ROI_profile = rowSums(is.finite(feature_matrix)) == length(feature_order),
  stringsAsFactors = FALSE
)
feature_data <- cbind(
  feature_data,
  as.data.frame(feature_matrix[participant_ids, feature_order, drop = FALSE], check.names = FALSE)
)
if ("calibration_status" %in% names(scores)) {
  feature_data$calibration_status <- vapply(
    participant_ids,
    function(id) collapse_unique(scores$calibration_status[scores$match_id == id]),
    character(1)
  )
  feature_data$contains_uncalibrated_site_score <- grepl(
    "uncalibrated", feature_data$calibration_status, ignore.case = TRUE
  )
}

# Optional CMD input is used only for Part 1 cognitive associations.
# It is not added to the Part 2/3 SVR predictor set.
cmd_available <- !identical(cmd_arg, "-")
cmd_file <- "not supplied"
if (cmd_available) {
  cmd_file <- normalizePath(cmd_arg, winslash = "/", mustWork = TRUE)
  cmd_input <- read.csv(cmd_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_cmd_columns <- c("match_id", "CMD")
  missing_cmd_columns <- setdiff(required_cmd_columns, names(cmd_input))
  if (length(missing_cmd_columns) > 0) {
    stop(
      "CMD input must be individual_CMD.csv from 05_CMD.R; missing: ",
      paste(missing_cmd_columns, collapse = ", ")
    )
  }
  cmd_input$match_id <- safe_chr(cmd_input$match_id)
  cmd_input$CMD_24ROI <- safe_num(cmd_input$CMD)
  cmd_input <- cmd_input[cmd_input$match_id != "", c("match_id", "CMD_24ROI"), drop = FALSE]
  if (anyDuplicated(cmd_input$match_id)) stop("CMD input contains duplicate match_id values")
  if (!any(is.finite(cmd_input$CMD_24ROI))) stop("CMD input contains no finite CMD values")

  feature_data$.feature_order <- seq_len(nrow(feature_data))
  feature_data <- merge(
    feature_data, cmd_input, by = "match_id", all.x = TRUE, sort = FALSE
  )
  feature_data <- feature_data[order(feature_data$.feature_order), , drop = FALSE]
  feature_data$.feature_order <- NULL
}

# ------------------------------------------------------------
# Prepare cognitive records and merge participant features.
# ------------------------------------------------------------
cognitive <- read_table_auto(cognitive_file)
cognitive$match_id <- make_match_id(cognitive)
required_cognitive_columns <- c("Age", "Sex", diagnosis_column)
missing_cognitive_columns <- setdiff(required_cognitive_columns, names(cognitive))
if (length(missing_cognitive_columns) > 0) {
  stop("Cognitive data are missing: ", paste(missing_cognitive_columns, collapse = ", "))
}
cognitive$Age <- safe_num(cognitive$Age)
cognitive$Sex <- safe_chr(cognitive$Sex)
cognitive$Diagnosis <- normalize_diagnosis(cognitive[[diagnosis_column]])
cognitive_outcomes <- trimws(unlist(strsplit(outcome_arg, ",", fixed = TRUE)))
cognitive_outcomes <- unique(cognitive_outcomes[cognitive_outcomes != ""])
available_outcomes <- intersect(cognitive_outcomes, names(cognitive))
if (length(available_outcomes) == 0) {
  stop("None of the requested cognitive outcome columns are present")
}
for (outcome in available_outcomes) cognitive[[outcome]] <- safe_num(cognitive[[outcome]])
has_follow_up_time <- follow_up_time_column %in% names(cognitive)
if (has_follow_up_time) {
  cognitive$follow_up_time_analysis <- safe_num(cognitive[[follow_up_time_column]])
} else {
  cognitive$follow_up_time_analysis <- NA_real_
}
cognitive$cognitive_row_id <- seq_len(nrow(cognitive))
cognitive <- cognitive[cognitive$match_id != "", , drop = FALSE]

analysis_records <- merge(cognitive, feature_data, by = "match_id", all.x = TRUE, sort = FALSE)
analysis_records <- analysis_records[
  analysis_records$complete_24ROI_profile %in% TRUE &
    is.finite(analysis_records$Age) & analysis_records$Sex != "" &
    analysis_records$Diagnosis != "",
  , drop = FALSE
]

healthy_labels <- normalize_diagnosis(trimws(unlist(strsplit(healthy_arg, ",", fixed = TRUE))))
healthy_labels <- healthy_labels[healthy_labels != ""]
if (length(healthy_labels) == 0) stop("At least one healthy-control label is required")
if (tolower(prediction_group_arg) == "auto") {
  prediction_groups <- sort(setdiff(unique(analysis_records$Diagnosis), healthy_labels))
} else {
  prediction_groups <- normalize_diagnosis(
    trimws(unlist(strsplit(prediction_group_arg, ",", fixed = TRUE)))
  )
  prediction_groups <- unique(prediction_groups[prediction_groups != ""])
}

task_log <- list()

select_cross_sectional_rows <- function(data, outcome) {
  candidates <- data[is.finite(data[[outcome]]), , drop = FALSE]
  if (has_follow_up_time) {
    baseline <- candidates[
      !is.finite(candidates$follow_up_time_analysis) |
        candidates$follow_up_time_analysis <= 0,
      , drop = FALSE
    ]
    if (nrow(baseline) > 0) candidates <- baseline
  }
  time_order <- ifelse(
    is.finite(candidates$follow_up_time_analysis),
    candidates$follow_up_time_analysis,
    0
  )
  candidates <- candidates[order(candidates$match_id, time_order, candidates$cognitive_row_id), , drop = FALSE]
  candidates[!duplicated(candidates$match_id), , drop = FALSE]
}

# ------------------------------------------------------------
# Part 1. Cognitive associations
# Partial Spearman correlation adjusted for age and sex.
# ------------------------------------------------------------
partial_spearman <- function(data, outcome, feature) {
  keep <- is.finite(data[[outcome]]) & is.finite(data[[feature]]) &
    is.finite(data$Age) & data$Sex != ""
  d <- data[keep, , drop = FALSE]
  base <- data.frame(
    n = nrow(d), partial_rho = NA_real_, p_value = NA_real_,
    status = "insufficient_data", stringsAsFactors = FALSE
  )
  if (nrow(d) < minimum_association_n ||
      stats::sd(d[[outcome]]) == 0 || stats::sd(d[[feature]]) == 0) {
    return(base)
  }
  sex_factor <- factor(d$Sex)
  if (nlevels(sex_factor) > 1) {
    sex_matrix <- stats::model.matrix(~ sex_factor)[, -1, drop = FALSE]
    covariate_matrix <- cbind(Intercept = 1, Age = d$Age, sex_matrix)
  } else {
    covariate_matrix <- cbind(Intercept = 1, Age = d$Age)
  }
  ranked_outcome <- rank(d[[outcome]], ties.method = "average")
  ranked_feature <- rank(d[[feature]], ties.method = "average")
  outcome_residual <- stats::lm.fit(covariate_matrix, ranked_outcome)$residuals
  feature_residual <- stats::lm.fit(covariate_matrix, ranked_feature)$residuals
  rho <- suppressWarnings(stats::cor(outcome_residual, feature_residual, method = "pearson"))
  degrees_freedom <- nrow(d) - qr(covariate_matrix)$rank - 1
  if (!is.finite(rho) || degrees_freedom < 1) return(base)
  rho_for_test <- max(min(rho, 1 - 1e-12), -1 + 1e-12)
  statistic <- rho_for_test * sqrt(degrees_freedom / (1 - rho_for_test^2))
  data.frame(
    n = nrow(d),
    partial_rho = rho,
    p_value = 2 * stats::pt(-abs(statistic), df = degrees_freedom),
    status = "completed",
    stringsAsFactors = FALSE
  )
}

association_rows <- list()
for (outcome in available_outcomes) {
  baseline_outcome <- select_cross_sectional_rows(analysis_records, outcome)
  association_groups <- list(
    All_participants = baseline_outcome,
    All_disorders_excluding_HC = baseline_outcome[
      !(baseline_outcome$Diagnosis %in% healthy_labels), , drop = FALSE
    ]
  )
  for (diagnosis in sort(unique(baseline_outcome$Diagnosis))) {
    association_groups[[paste0("Diagnosis_", diagnosis)]] <- baseline_outcome[
      baseline_outcome$Diagnosis == diagnosis, , drop = FALSE
    ]
  }
  for (group_name in names(association_groups)) {
    group_data <- association_groups[[group_name]]
    for (feature in feature_order) {
      result <- partial_spearman(group_data, outcome, feature)
      association_rows[[length(association_rows) + 1]] <- data.frame(
        analysis_part = "Part_1_Cognitive_associations",
        analysis_group = group_name,
        cognitive_outcome = outcome,
        feature_family = "regional_centile",
        feature = feature,
        result,
        method = "partial Spearman correlation adjusted for age and sex",
        stringsAsFactors = FALSE
      )
    }
    if (cmd_available) {
      result <- partial_spearman(group_data, outcome, "CMD_24ROI")
      association_rows[[length(association_rows) + 1]] <- data.frame(
        analysis_part = "Part_1_Cognitive_associations",
        analysis_group = group_name,
        cognitive_outcome = outcome,
        feature_family = "CMD",
        feature = "CMD_24ROI",
        result,
        method = "partial Spearman correlation adjusted for age and sex",
        stringsAsFactors = FALSE
      )
    }
  }
}
associations <- bind_rows_base(association_rows)
associations$pFDR <- NA_real_
if (nrow(associations) > 0) {
  regional_index <- associations$feature_family == "regional_centile"
  fdr_groups <- unique(paste(
    associations$analysis_group[regional_index],
    associations$cognitive_outcome[regional_index],
    sep = "\r"
  ))
  for (fdr_group in fdr_groups) {
    index <- which(regional_index & paste(
      associations$analysis_group, associations$cognitive_outcome, sep = "\r"
    ) == fdr_group)
    associations$pFDR[index] <- stats::p.adjust(associations$p_value[index], method = "BH")
  }

  cmd_index <- associations$feature_family == "CMD"
  if (any(cmd_index)) {
    for (group_name in unique(associations$analysis_group[cmd_index])) {
      index <- which(cmd_index & associations$analysis_group == group_name)
      associations$pFDR[index] <- stats::p.adjust(
        associations$p_value[index], method = "BH"
      )
    }
  }

  write.csv(
    associations,
    file.path(output_dir, "Part1_cognitive_associations_partial_Spearman.csv"),
    row.names = FALSE
  )
  if (any(cmd_index)) {
    write.csv(
      associations[cmd_index, , drop = FALSE],
      file.path(output_dir, "Part1_CMD_cognitive_associations_partial_Spearman.csv"),
      row.names = FALSE
    )
  }
}

# ------------------------------------------------------------
# Shared explicit nested-CV SVR engine for Parts 2 and 3.
# ------------------------------------------------------------
make_predictor_matrix <- function(data) {
  numeric_predictors <- data.frame(Age = data$Age, data[, feature_order, drop = FALSE], check.names = FALSE)
  sex_factor <- factor(data$Sex)
  if (nlevels(sex_factor) > 1) {
    sex_matrix <- stats::model.matrix(~ sex_factor)[, -1, drop = FALSE]
    colnames(sex_matrix) <- paste0("Sex_", make.names(levels(sex_factor)[-1]))
    numeric_predictors <- cbind(numeric_predictors, sex_matrix)
  } else {
    numeric_predictors$Sex_single_level <- 0
  }
  as.data.frame(numeric_predictors, check.names = FALSE)
}

make_outer_test_folds <- function(data, outcome, group_by_subject) {
  if (!group_by_subject) {
    set.seed(random_seed)
    return(caret::createFolds(data[[outcome]], k = outer_folds, returnTrain = FALSE))
  }
  subject_outcome <- tapply(data[[outcome]], data$match_id, mean, na.rm = TRUE)
  set.seed(random_seed)
  subject_folds <- caret::createFolds(
    as.numeric(subject_outcome), k = outer_folds, returnTrain = FALSE
  )
  lapply(subject_folds, function(subject_index) {
    which(data$match_id %in% names(subject_outcome)[subject_index])
  })
}

prediction_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  if (length(observed) < 4 || stats::sd(observed) == 0 || stats::sd(predicted) == 0) {
    return(data.frame(
      n = length(observed), Pearson_r = NA_real_, Pearson_r_95CI_low = NA_real_,
      Pearson_r_95CI_high = NA_real_, p_value = NA_real_, RMSE = NA_real_,
      MAE = NA_real_, R2 = NA_real_, stringsAsFactors = FALSE
    ))
  }
  correlation <- suppressWarnings(stats::cor.test(observed, predicted, method = "pearson"))
  data.frame(
    n = length(observed),
    Pearson_r = unname(correlation$estimate),
    Pearson_r_95CI_low = correlation$conf.int[[1]],
    Pearson_r_95CI_high = correlation$conf.int[[2]],
    p_value = correlation$p.value,
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted)),
    R2 = 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2),
    stringsAsFactors = FALSE
  )
}

run_nested_svr <- function(data, outcome, analysis_group, analysis_part, group_by_subject) {
  data <- data[
    is.finite(data[[outcome]]) & is.finite(data$Age) & data$Sex != "" &
      rowSums(is.finite(as.matrix(data[, feature_order, drop = FALSE]))) == length(feature_order),
    , drop = FALSE
  ]
  if (!group_by_subject) {
    data <- data[order(data$match_id, data$cognitive_row_id), , drop = FALSE]
    data <- data[!duplicated(data$match_id), , drop = FALSE]
  }
  n_subjects <- length(unique(data$match_id))
  if (nrow(data) < minimum_prediction_n || n_subjects < outer_folds) {
    return(list(
      predictions = data.frame(), tuning = data.frame(), metrics = data.frame(),
      status = "skipped", reason = "insufficient observations or subjects for 10-fold CV"
    ))
  }

  predictor_matrix <- make_predictor_matrix(data)
  observed <- data[[outcome]]
  test_folds <- make_outer_test_folds(data, outcome, group_by_subject)
  prediction_rows <- list()
  tuning_rows <- list()

  for (fold_index in seq_along(test_folds)) {
    test_index <- test_folds[[fold_index]]
    training_index <- setdiff(seq_len(nrow(data)), test_index)
    x_train <- predictor_matrix[training_index, , drop = FALSE]
    x_test <- predictor_matrix[test_index, , drop = FALSE]
    y_train <- observed[training_index]

    zero_variance <- caret::nearZeroVar(x_train, saveMetrics = TRUE)
    keep_predictors <- rownames(zero_variance)[
      !zero_variance$zeroVar & !zero_variance$nzv
    ]
    if (length(keep_predictors) < 1) {
      return(list(
        predictions = data.frame(), tuning = data.frame(), metrics = data.frame(),
        status = "failed", reason = paste0("no usable predictors in outer fold ", fold_index)
      ))
    }
    x_train <- x_train[, keep_predictors, drop = FALSE]
    x_test <- x_test[, keep_predictors, drop = FALSE]

    set.seed(random_seed + fold_index)
    inner_control <- caret::trainControl(
      method = "cv",
      number = inner_folds,
      allowParallel = FALSE
    )
    fit_error <- NULL
    fit <- tryCatch(
      caret::train(
        x = x_train,
        y = y_train,
        method = "svmRadial",
        metric = "RMSE",
        trControl = inner_control,
        tuneLength = 10,
        preProcess = c("center", "scale")
      ),
      error = function(error) {
        fit_error <<- conditionMessage(error)
        NULL
      }
    )
    if (is.null(fit)) {
      return(list(
        predictions = data.frame(), tuning = data.frame(), metrics = data.frame(),
        status = "failed",
        reason = paste0("SVR failed in outer fold ", fold_index, ": ", fit_error)
      ))
    }
    fold_prediction <- as.numeric(stats::predict(fit, newdata = x_test))
    prediction_rows[[length(prediction_rows) + 1]] <- data.frame(
      analysis_part = analysis_part,
      analysis_group = analysis_group,
      cognitive_outcome = outcome,
      match_id = data$match_id[test_index],
      cognitive_row_id = data$cognitive_row_id[test_index],
      Diagnosis = data$Diagnosis[test_index],
      outer_fold = fold_index,
      observed = observed[test_index],
      predicted = fold_prediction,
      stringsAsFactors = FALSE
    )
    tuning_rows[[length(tuning_rows) + 1]] <- data.frame(
      analysis_part = analysis_part,
      analysis_group = analysis_group,
      cognitive_outcome = outcome,
      outer_fold = fold_index,
      training_n = length(training_index),
      heldout_n = length(test_index),
      n_predictors_used = length(keep_predictors),
      best_sigma = fit$bestTune$sigma[[1]],
      best_C = fit$bestTune$C[[1]],
      stringsAsFactors = FALSE
    )
  }

  predictions <- bind_rows_base(prediction_rows)
  metrics <- data.frame(
    analysis_part = analysis_part,
    analysis_group = analysis_group,
    cognitive_outcome = outcome,
    prediction_metrics(predictions$observed, predictions$predicted),
    outer_CV = "10-fold",
    tuning = "5-fold CV within each outer training set",
    model = "radial-basis-function support vector regression",
    stringsAsFactors = FALSE
  )
  list(
    predictions = predictions,
    tuning = bind_rows_base(tuning_rows),
    metrics = metrics,
    status = "completed",
    reason = ""
  )
}

# ------------------------------------------------------------
# Part 2. Cross-sectional cognitive prediction
# ------------------------------------------------------------
cross_prediction_rows <- list()
cross_tuning_rows <- list()
cross_metric_rows <- list()

for (diagnosis in prediction_groups) {
  for (outcome in available_outcomes) {
    baseline_outcome <- select_cross_sectional_rows(analysis_records, outcome)
    task_data <- baseline_outcome[baseline_outcome$Diagnosis == diagnosis, , drop = FALSE]
    task_name <- paste(diagnosis, outcome, sep = "__")
    message("Cross-sectional SVR: ", task_name)
    result <- run_nested_svr(
      task_data, outcome, diagnosis,
      "Part_2_Cross_sectional_cognitive_prediction", FALSE
    )
    task_log[[length(task_log) + 1]] <- data.frame(
      analysis_part = "Part_2_Cross_sectional_cognitive_prediction",
      analysis_group = diagnosis,
      cognitive_outcome = outcome,
      status = result$status,
      reason = result$reason,
      stringsAsFactors = FALSE
    )
    cross_prediction_rows[[length(cross_prediction_rows) + 1]] <- result$predictions
    cross_tuning_rows[[length(cross_tuning_rows) + 1]] <- result$tuning
    cross_metric_rows[[length(cross_metric_rows) + 1]] <- result$metrics
  }
}

cross_predictions <- bind_rows_base(cross_prediction_rows)
cross_tuning <- bind_rows_base(cross_tuning_rows)
cross_metrics <- bind_rows_base(cross_metric_rows)
if (nrow(cross_metrics) > 0) {
  cross_metrics$pFDR <- stats::p.adjust(cross_metrics$p_value, method = "BH")
  write.csv(
    cross_metrics,
    file.path(output_dir, "Part2_cross_sectional_SVR_10fold_metrics.csv"),
    row.names = FALSE
  )
  write.csv(
    cross_predictions,
    file.path(output_dir, "Part2_cross_sectional_SVR_10fold_heldout_predictions.csv"),
    row.names = FALSE
  )
  write.csv(
    cross_tuning,
    file.path(output_dir, "Part2_cross_sectional_SVR_outer_fold_tuning.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# Part 3. Longitudinal cognitive prediction
# ------------------------------------------------------------
long_prediction_rows <- list()
long_tuning_rows <- list()
long_metric_rows <- list()

if (has_follow_up_time) {
  longitudinal_records <- analysis_records[
    is.finite(analysis_records$follow_up_time_analysis) &
      analysis_records$follow_up_time_analysis > 0,
    , drop = FALSE
  ]
  for (diagnosis in prediction_groups) {
    for (outcome in available_outcomes) {
      task_data <- longitudinal_records[
        longitudinal_records$Diagnosis == diagnosis &
          is.finite(longitudinal_records[[outcome]]),
        , drop = FALSE
      ]
      task_name <- paste(diagnosis, outcome, sep = "__")
      message("Longitudinal SVR: ", task_name)
      result <- run_nested_svr(
        task_data, outcome, diagnosis,
        "Part_3_Longitudinal_cognitive_prediction", TRUE
      )
      task_log[[length(task_log) + 1]] <- data.frame(
        analysis_part = "Part_3_Longitudinal_cognitive_prediction",
        analysis_group = diagnosis,
        cognitive_outcome = outcome,
        status = result$status,
        reason = result$reason,
        stringsAsFactors = FALSE
      )
      long_prediction_rows[[length(long_prediction_rows) + 1]] <- result$predictions
      long_tuning_rows[[length(long_tuning_rows) + 1]] <- result$tuning
      long_metric_rows[[length(long_metric_rows) + 1]] <- result$metrics
    }
  }
} else {
  task_log[[length(task_log) + 1]] <- data.frame(
    analysis_part = "Part_3_Longitudinal_cognitive_prediction",
    analysis_group = "All",
    cognitive_outcome = "All",
    status = "skipped",
    reason = paste0("follow-up time column not found: ", follow_up_time_column),
    stringsAsFactors = FALSE
  )
}

long_predictions <- bind_rows_base(long_prediction_rows)
long_tuning <- bind_rows_base(long_tuning_rows)
long_metrics <- bind_rows_base(long_metric_rows)
if (nrow(long_metrics) > 0) {
  long_metrics$pFDR <- stats::p.adjust(long_metrics$p_value, method = "BH")
  write.csv(
    long_metrics,
    file.path(output_dir, "Part3_longitudinal_SVR_10fold_metrics.csv"),
    row.names = FALSE
  )
  write.csv(
    long_predictions,
    file.path(output_dir, "Part3_longitudinal_SVR_10fold_heldout_predictions.csv"),
    row.names = FALSE
  )
  write.csv(
    long_tuning,
    file.path(output_dir, "Part3_longitudinal_SVR_outer_fold_tuning.csv"),
    row.names = FALSE
  )
}

task_log_table <- bind_rows_base(task_log)
if (nrow(task_log_table) > 0) {
  write.csv(
    task_log_table,
    file.path(output_dir, "cognitive_analysis_task_log.csv"),
    row.names = FALSE
  )
}

settings <- data.frame(
  item = c(
    "score input", "CMD input", "cognitive input", "cognitive outcomes requested",
    "cognitive outcomes available", "follow-up time column", "feature set",
    "Part 1 method", "regional association FDR family", "CMD association FDR family",
    "prediction model", "outer CV", "inner tuning",
    "standardization", "prediction metric", "longitudinal fold grouping",
    "minimum association n", "minimum prediction n", "seed",
    "input percentile scale"
  ),
  value = c(
    score_file,
    cmd_file,
    cognitive_file,
    paste(cognitive_outcomes, collapse = ","),
    paste(available_outcomes, collapse = ","),
    if (has_follow_up_time) follow_up_time_column else "not available",
    "Age, Sex, and 24 bilateral hippocampal percentile features",
    "partial Spearman correlation adjusted for age and sex",
    "within each analysis group and cognitive outcome across 24 regional centiles",
    if (cmd_available) {
      "within each analysis group across available cognitive outcomes"
    } else {
      "not applicable; CMD input not supplied"
    },
    "radial-basis-function support vector regression",
    "10-fold held-out predictions",
    "5-fold CV within each outer training set; tuneLength=10",
    "centering and scaling fitted using training data only",
    "Pearson correlation across pooled held-out predictions",
    "all rows from one match_id remain in the same outer fold",
    minimum_association_n,
    minimum_prediction_n,
    random_seed,
    percentile_scale
  ),
  stringsAsFactors = FALSE
)
write.csv(settings, file.path(output_dir, "cognitive_analysis_settings.csv"), row.names = FALSE)

message(
  "Cognitive analyses completed. Associations: ", nrow(associations),
  "; CMD-cognition associations: ", sum(associations$feature_family == "CMD"),
  "; cross-sectional prediction tasks: ", nrow(cross_metrics),
  "; longitudinal prediction tasks: ", nrow(long_metrics),
  ". Outputs: ", output_dir
)

# ============================================================
# 06_Classification_DPS.R
# Held-out disease classification and disease propensity scores
# from 24 bilateral hippocampal normative percentile features.
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Usage:
# Rscript 06_Classification_DPS.R <normative_scores_long.csv> \
#   <clinical_data.csv> <output_dir> [disease_labels|auto] \
#   [minimum_group_n] [seed] [diagnosis_column] [healthy_control_labels]
#
# Default disease labels: MCI,AD,PD,CSVD,MS
# Comma-separate multiple disease or healthy-control labels.
#
# Method implemented here:
#   1. Each HC-versus-disease dataset is randomly split 70:30 into
#      training and independent held-out test partitions.
#   2. Within each partition separately, disease cases are matched
#      1:1 to same-sex HCs by nearest age, without replacement.
#   3. Regularized logistic regression and a single-hidden-layer
#      neural network are tuned by 5-fold CV in the training set.
#   4. The candidate with the highest training CV ROC AUC is selected.
#   5. AUC and DPS are calculated only from held-out probabilities.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 06_Classification_DPS.R <normative_scores_long.csv> ",
    "<clinical_data.csv> <output_dir> [disease_labels|auto] ",
    "[minimum_group_n] [seed] [diagnosis_column] [healthy_control_labels]"
  )
}

score_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
disease_arg <- if (length(args) >= 4) args[[4]] else "MCI,AD,PD,CSVD,MS"
minimum_group_n <- if (length(args) >= 5) as.integer(args[[5]]) else 50L
random_seed <- if (length(args) >= 6) as.integer(args[[6]]) else 20260520L
diagnosis_column <- if (length(args) >= 7) args[[7]] else "Diagnosis"
healthy_arg <- if (length(args) >= 8) args[[8]] else "HC"

if (!is.finite(minimum_group_n) || minimum_group_n < 10) {
  stop("minimum_group_n must be an integer of at least 10")
}
if (!is.finite(random_seed)) stop("seed must be an integer")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

required_packages <- c("caret", "pROC", "glmnet", "nnet")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

set.seed(random_seed)
training_fraction <- 0.70
test_fraction <- 0.30
cv_folds <- 5L

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
    "Clinical data require either match_id or the same FreeSurfer path columns ",
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

roi_order <- c(
  "Whole_hippocampus", "subiculum", "CA1", "presubiculum",
  "molecular_layer_HP", "GC.ML.DG", "CA3", "CA4",
  "Hippocampal_tail", "parasubiculum", "fimbria", "HATA"
)
feature_order <- unlist(
  lapply(c("left", "right"), function(side) paste0(side, "_", roi_order)),
  use.names = FALSE
)

clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
clinical$match_id <- make_match_id(clinical)
missing_clinical <- setdiff(c("Age", "Sex", diagnosis_column), names(clinical))
if (length(missing_clinical) > 0) {
  stop("Clinical data are missing: ", paste(missing_clinical, collapse = ", "))
}
clinical$Age <- safe_num(clinical$Age)
clinical$Sex <- safe_chr(clinical$Sex)
clinical$Diagnosis <- normalize_diagnosis(clinical[[diagnosis_column]])
clinical <- clinical[clinical$match_id != "", , drop = FALSE]
if (anyDuplicated(clinical$match_id)) stop("Clinical data contain duplicate match_id values")
clinical <- clinical[, c("match_id", "Age", "Sex", "Diagnosis"), drop = FALSE]

healthy_labels <- normalize_diagnosis(trimws(unlist(strsplit(healthy_arg, ",", fixed = TRUE))))
healthy_labels <- healthy_labels[healthy_labels != ""]
if (length(healthy_labels) == 0) stop("At least one healthy-control label is required")

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
scores$hemisphere_classification <- normalize_hemisphere(scores$hemisphere)
scores$ROI_classification <- normalize_roi(scores$ROI)
if (any(scores$match_id == "")) stop("The score file contains empty match_id values")
if (anyNA(scores$hemisphere_classification)) {
  stop("Unrecognized hemisphere values in score input")
}
if (anyNA(scores$ROI_classification)) stop("Unrecognized ROI values in score input")
scores$feature_name_classification <- paste0(
  scores$hemisphere_classification, "_", scores$ROI_classification
)
unexpected_features <- setdiff(unique(scores$feature_name_classification), feature_order)
if (length(unexpected_features) > 0) {
  stop("Unexpected classification features: ", paste(unexpected_features, collapse = ", "))
}
duplicate_key <- duplicated(scores[c("match_id", "feature_name_classification")]) |
  duplicated(scores[c("match_id", "feature_name_classification")], fromLast = TRUE)
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
  match(scores$feature_name_classification, feature_order)
)] <- percentile

participant_audit <- data.frame(
  match_id = participant_ids,
  n_hippocampal_features = rowSums(is.finite(feature_matrix)),
  complete_24ROI_profile = rowSums(is.finite(feature_matrix)) == length(feature_order),
  stringsAsFactors = FALSE
)
if ("calibration_status" %in% names(scores)) {
  participant_audit$calibration_status <- vapply(
    participant_ids,
    function(id) collapse_unique(scores$calibration_status[scores$match_id == id]),
    character(1)
  )
  participant_audit$contains_uncalibrated_site_score <- grepl(
    "uncalibrated", participant_audit$calibration_status, ignore.case = TRUE
  )
}

analysis_data <- merge(participant_audit, clinical, by = "match_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(participant_ids, analysis_data$match_id), , drop = FALSE]
analysis_data <- cbind(
  analysis_data,
  as.data.frame(feature_matrix[analysis_data$match_id, feature_order, drop = FALSE], check.names = FALSE)
)
analysis_data$row_id <- seq_len(nrow(analysis_data))
analysis_data <- analysis_data[
  analysis_data$complete_24ROI_profile & is.finite(analysis_data$Age) &
    analysis_data$Sex != "" & analysis_data$Diagnosis != "",
  , drop = FALSE
]

if (tolower(disease_arg) == "auto") {
  disease_labels <- sort(setdiff(unique(analysis_data$Diagnosis), healthy_labels))
} else {
  disease_labels <- normalize_diagnosis(trimws(unlist(strsplit(disease_arg, ",", fixed = TRUE))))
  disease_labels <- unique(disease_labels[disease_labels != ""])
}
if (length(disease_labels) == 0) stop("No disease labels were requested or detected")

diagnosis_counts <- as.data.frame(table(analysis_data$Diagnosis), stringsAsFactors = FALSE)
names(diagnosis_counts) <- c("Diagnosis", "n_complete_24ROI")
eligible_diseases <- disease_labels[
  vapply(disease_labels, function(label) {
    sum(analysis_data$Diagnosis == label) >= minimum_group_n
  }, logical(1))
]

skips <- list()
for (label in setdiff(disease_labels, eligible_diseases)) {
  skips[[length(skips) + 1]] <- data.frame(
    disease = label,
    stage = "eligibility",
    reason = paste0("fewer than minimum_group_n = ", minimum_group_n, " complete cases"),
    stringsAsFactors = FALSE
  )
}
if (length(eligible_diseases) == 0) {
  stop("No disease group met minimum_group_n = ", minimum_group_n)
}
if (sum(analysis_data$Diagnosis %in% healthy_labels) < minimum_group_n) {
  stop("Healthy controls have fewer than minimum_group_n complete cases")
}

glmnet_grid <- expand.grid(
  alpha = c(0, 0.25, 0.5, 0.75, 1),
  lambda = 10^seq(-4, 1, length.out = 20)
)
nnet_grid <- expand.grid(
  size = c(1, 3, 5, 7),
  decay = c(0, 0.001, 0.01, 0.1)
)

match_partition <- function(partition, disease_label, task_name, partition_name) {
  cases <- partition[partition$Diagnosis == disease_label, , drop = FALSE]
  controls <- partition[partition$Diagnosis %in% healthy_labels, , drop = FALSE]
  cases <- cases[order(cases$Sex, cases$Age, cases$match_id), , drop = FALSE]
  available_controls <- controls
  matched_case_rows <- list()
  matched_control_rows <- list()
  detail_rows <- list()

  for (index in seq_len(nrow(cases))) {
    current_case <- cases[index, , drop = FALSE]
    candidates <- which(available_controls$Sex == current_case$Sex[[1]])
    if (length(candidates) == 0) next
    age_difference <- abs(available_controls$Age[candidates] - current_case$Age[[1]])
    candidate_order <- order(age_difference, available_controls$match_id[candidates])
    selected_position <- candidates[candidate_order[[1]]]
    selected_control <- available_controls[selected_position, , drop = FALSE]
    pair_id <- paste(
      task_name, partition_name, current_case$match_id[[1]],
      selected_control$match_id[[1]], sep = "__"
    )
    difference <- abs(selected_control$Age[[1]] - current_case$Age[[1]])
    current_case$matched_role <- "Case"
    current_case$match_pair_id <- pair_id
    current_case$match_age_difference <- difference
    selected_control$matched_role <- "Matched_HC"
    selected_control$match_pair_id <- pair_id
    selected_control$match_age_difference <- difference
    matched_case_rows[[length(matched_case_rows) + 1]] <- current_case
    matched_control_rows[[length(matched_control_rows) + 1]] <- selected_control
    detail_rows[[length(detail_rows) + 1]] <- data.frame(
      task = task_name,
      partition = partition_name,
      pair_id = pair_id,
      case_match_id = current_case$match_id[[1]],
      HC_match_id = selected_control$match_id[[1]],
      Sex = current_case$Sex[[1]],
      case_Age = current_case$Age[[1]],
      HC_Age = selected_control$Age[[1]],
      absolute_age_difference = difference,
      stringsAsFactors = FALSE
    )
    available_controls <- available_controls[-selected_position, , drop = FALSE]
  }

  matched_cases <- bind_rows_base(matched_case_rows)
  matched_controls <- bind_rows_base(matched_control_rows)
  matched <- bind_rows_base(list(matched_cases, matched_controls))
  if (nrow(matched) > 0) {
    matched$outcome <- factor(
      ifelse(matched$Diagnosis == disease_label, "Case", "Reference"),
      levels = c("Case", "Reference")
    )
  }
  list(
    matched = matched,
    details = bind_rows_base(detail_rows),
    n_cases_before_matching = nrow(cases),
    n_controls_before_matching = nrow(controls),
    n_matched_pairs = nrow(matched_cases),
    n_unmatched_cases = nrow(cases) - nrow(matched_cases)
  )
}

best_cv_row <- function(fit) {
  candidate <- fit$results
  for (parameter in names(fit$bestTune)) {
    candidate <- candidate[candidate[[parameter]] == fit$bestTune[[parameter]], , drop = FALSE]
  }
  candidate[which.max(candidate$ROC), , drop = FALSE]
}

fit_candidate <- function(x_train, y_train, model_name, train_control) {
  if (model_name == "regularized_logistic_regression") {
    caret::train(
      x = x_train,
      y = y_train,
      method = "glmnet",
      metric = "ROC",
      trControl = train_control,
      tuneGrid = glmnet_grid,
      preProcess = c("center", "scale"),
      family = "binomial"
    )
  } else if (model_name == "single_hidden_layer_neural_network") {
    caret::train(
      x = x_train,
      y = y_train,
      method = "nnet",
      metric = "ROC",
      trControl = train_control,
      tuneGrid = nnet_grid,
      preProcess = c("center", "scale"),
      trace = FALSE,
      MaxNWts = 5000,
      maxit = 300
    )
  } else {
    stop("Unknown candidate model: ", model_name)
  }
}

candidate_names <- c(
  "regularized_logistic_regression",
  "single_hidden_layer_neural_network"
)

metrics_list <- list()
cv_list <- list()
prediction_list <- list()
roc_list <- list()
confusion_list <- list()
matching_list <- list()
partition_list <- list()
selected_models <- list()

for (task_index in seq_along(eligible_diseases)) {
  disease <- eligible_diseases[[task_index]]
  task_name <- paste0("HC_vs_", disease)
  message("Classification task: ", task_name)
  task_data <- analysis_data[
    analysis_data$Diagnosis %in% c(healthy_labels, disease), , drop = FALSE
  ]
  task_data$outcome_for_split <- factor(
    ifelse(task_data$Diagnosis == disease, "Case", "Reference"),
    levels = c("Case", "Reference")
  )

  set.seed(random_seed + task_index - 1L)
  training_index <- as.vector(caret::createDataPartition(
    task_data$outcome_for_split, p = training_fraction, list = FALSE
  ))
  training_partition <- task_data[training_index, , drop = FALSE]
  test_partition <- task_data[-training_index, , drop = FALSE]

  training_match <- match_partition(training_partition, disease, task_name, "training")
  test_match <- match_partition(test_partition, disease, task_name, "heldout_test")
  matching_list[[length(matching_list) + 1]] <- training_match$details
  matching_list[[length(matching_list) + 1]] <- test_match$details
  partition_list[[length(partition_list) + 1]] <- data.frame(
    task = task_name,
    partition = c("training", "heldout_test"),
    n_cases_before_matching = c(
      training_match$n_cases_before_matching, test_match$n_cases_before_matching
    ),
    n_HC_before_matching = c(
      training_match$n_controls_before_matching, test_match$n_controls_before_matching
    ),
    n_matched_pairs = c(training_match$n_matched_pairs, test_match$n_matched_pairs),
    n_unmatched_cases = c(training_match$n_unmatched_cases, test_match$n_unmatched_cases),
    stringsAsFactors = FALSE
  )

  training_data <- training_match$matched
  test_data <- test_match$matched
  if (nrow(training_data) == 0 || nrow(test_data) == 0 ||
      min(table(training_data$outcome)) < cv_folds || min(table(test_data$outcome)) < 2) {
    skips[[length(skips) + 1]] <- data.frame(
      disease = disease,
      stage = "partition_matching",
      reason = "insufficient matched observations for 5-fold CV or held-out AUC",
      stringsAsFactors = FALSE
    )
    next
  }

  x_train <- training_data[, feature_order, drop = FALSE]
  x_test <- test_data[, feature_order, drop = FALSE]
  y_train <- training_data$outcome
  train_control <- caret::trainControl(
    method = "cv",
    number = cv_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )

  model_fits <- list()
  cv_task_rows <- list()
  for (model_index in seq_along(candidate_names)) {
    model_name <- candidate_names[[model_index]]
    set.seed(random_seed + task_index * 100L + model_index)
    fit_error <- NULL
    fit <- tryCatch(
      fit_candidate(x_train, y_train, model_name, train_control),
      error = function(error) {
        fit_error <<- conditionMessage(error)
        NULL
      }
    )
    if (is.null(fit)) {
      cv_task_rows[[length(cv_task_rows) + 1]] <- data.frame(
        task = task_name,
        disease = disease,
        candidate_model = model_name,
        cv_ROC_AUC = NA_real_,
        cv_sensitivity = NA_real_,
        cv_specificity = NA_real_,
        best_tuning_parameters = "",
        fit_status = "failed",
        error_message = fit_error,
        stringsAsFactors = FALSE
      )
      next
    }
    model_fits[[model_name]] <- fit
    best <- best_cv_row(fit)
    cv_task_rows[[length(cv_task_rows) + 1]] <- data.frame(
      task = task_name,
      disease = disease,
      candidate_model = model_name,
      cv_ROC_AUC = best$ROC[[1]],
      cv_sensitivity = best$Sens[[1]],
      cv_specificity = best$Spec[[1]],
      best_tuning_parameters = paste(
        names(fit$bestTune), unlist(fit$bestTune), sep = "=", collapse = ";"
      ),
      fit_status = "completed",
      error_message = "",
      stringsAsFactors = FALSE
    )
  }
  cv_task <- bind_rows_base(cv_task_rows)
  cv_list[[length(cv_list) + 1]] <- cv_task
  eligible_models <- cv_task[
    cv_task$fit_status == "completed" & is.finite(cv_task$cv_ROC_AUC), , drop = FALSE
  ]
  if (nrow(eligible_models) == 0) {
    skips[[length(skips) + 1]] <- data.frame(
      disease = disease,
      stage = "candidate_model_fitting",
      reason = "all candidate models failed",
      stringsAsFactors = FALSE
    )
    next
  }

  selected_name <- eligible_models$candidate_model[which.max(eligible_models$cv_ROC_AUC)]
  selected_model <- model_fits[[selected_name]]
  heldout_probability <- as.numeric(
    predict(selected_model, newdata = x_test, type = "prob")[, "Case"]
  )
  heldout_prediction <- factor(
    ifelse(heldout_probability >= 0.5, "Case", "Reference"),
    levels = c("Case", "Reference")
  )
  roc_object <- pROC::roc(
    response = test_data$outcome,
    predictor = heldout_probability,
    levels = c("Reference", "Case"),
    direction = "<",
    quiet = TRUE
  )
  auc_ci <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))
  confusion <- caret::confusionMatrix(
    heldout_prediction, test_data$outcome, positive = "Case"
  )

  selected_models[[task_name]] <- list(
    disease = disease,
    selected_model_name = selected_name,
    feature_order = feature_order,
    fitted_model = selected_model
  )
  metrics_list[[length(metrics_list) + 1]] <- data.frame(
    task = task_name,
    disease = disease,
    selected_model = selected_name,
    training_n = nrow(training_data),
    training_case_n = sum(training_data$outcome == "Case"),
    training_HC_n = sum(training_data$outcome == "Reference"),
    heldout_test_n = nrow(test_data),
    heldout_case_n = sum(test_data$outcome == "Case"),
    heldout_HC_n = sum(test_data$outcome == "Reference"),
    selected_training_CV_AUC = max(eligible_models$cv_ROC_AUC),
    heldout_AUC = as.numeric(pROC::auc(roc_object)),
    heldout_AUC_95CI_low = auc_ci[[1]],
    heldout_AUC_95CI_mid = auc_ci[[2]],
    heldout_AUC_95CI_high = auc_ci[[3]],
    heldout_accuracy = unname(confusion$overall["Accuracy"]),
    heldout_balanced_accuracy = unname(confusion$byClass["Balanced Accuracy"]),
    heldout_sensitivity = unname(confusion$byClass["Sensitivity"]),
    heldout_specificity = unname(confusion$byClass["Specificity"]),
    heldout_PPV = unname(confusion$byClass["Pos Pred Value"]),
    heldout_NPV = unname(confusion$byClass["Neg Pred Value"]),
    stringsAsFactors = FALSE
  )

  prediction_columns <- c(
    "match_id", "Diagnosis", "Age", "Sex", "matched_role",
    "match_pair_id", "match_age_difference"
  )
  if ("calibration_status" %in% names(test_data)) {
    prediction_columns <- c(prediction_columns, "calibration_status")
  }
  if ("contains_uncalibrated_site_score" %in% names(test_data)) {
    prediction_columns <- c(prediction_columns, "contains_uncalibrated_site_score")
  }
  prediction_output <- test_data[, prediction_columns, drop = FALSE]
  prediction_output$task <- task_name
  prediction_output$target_disease <- disease
  prediction_output$selected_model <- selected_name
  prediction_output$truth <- as.character(test_data$outcome)
  prediction_output$predicted_class <- as.character(heldout_prediction)
  prediction_output$DPS <- heldout_probability
  prediction_list[[length(prediction_list) + 1]] <- prediction_output

  roc_list[[length(roc_list) + 1]] <- data.frame(
    task = task_name,
    disease = disease,
    selected_model = selected_name,
    false_positive_rate = 1 - roc_object$specificities,
    sensitivity = roc_object$sensitivities,
    stringsAsFactors = FALSE
  )
  confusion_table <- as.data.frame(confusion$table, stringsAsFactors = FALSE)
  names(confusion_table) <- c("predicted", "truth", "n")
  confusion_table$task <- task_name
  confusion_table$disease <- disease
  confusion_list[[length(confusion_list) + 1]] <- confusion_table
}

metrics <- bind_rows_base(metrics_list)
cv_results <- bind_rows_base(cv_list)
predictions <- bind_rows_base(prediction_list)
roc_coordinates <- bind_rows_base(roc_list)
confusion_matrices <- bind_rows_base(confusion_list)
matching_details <- bind_rows_base(matching_list)
partition_audit <- bind_rows_base(partition_list)
skip_log <- bind_rows_base(skips)

write.csv(diagnosis_counts, file.path(output_dir, "classification_diagnosis_counts.csv"), row.names = FALSE)
if (nrow(metrics) > 0) {
  write.csv(metrics, file.path(output_dir, "classification_heldout_metrics.csv"), row.names = FALSE)
}
if (nrow(cv_results) > 0) {
  write.csv(
    cv_results,
    file.path(output_dir, "classification_training_CV_model_comparison.csv"),
    row.names = FALSE
  )
}
if (nrow(predictions) > 0) {
  write.csv(
    predictions,
    file.path(output_dir, "classification_heldout_predictions_DPS.csv"),
    row.names = FALSE
  )
}
if (nrow(roc_coordinates) > 0) {
  write.csv(
    roc_coordinates,
    file.path(output_dir, "classification_heldout_ROC_coordinates.csv"),
    row.names = FALSE
  )
}
if (nrow(confusion_matrices) > 0) {
  write.csv(
    confusion_matrices,
    file.path(output_dir, "classification_heldout_confusion_matrices.csv"),
    row.names = FALSE
  )
}
if (nrow(matching_details) > 0) {
  write.csv(
    matching_details,
    file.path(output_dir, "classification_age_sex_matching_details.csv"),
    row.names = FALSE
  )
}
if (nrow(partition_audit) > 0) {
  write.csv(
    partition_audit,
    file.path(output_dir, "classification_partition_matching_audit.csv"),
    row.names = FALSE
  )
}
if (nrow(skip_log) > 0) {
  write.csv(skip_log, file.path(output_dir, "classification_skip_log.csv"), row.names = FALSE)
}
if (length(selected_models) > 0) {
  saveRDS(selected_models, file.path(output_dir, "classification_selected_models.rds"))
}

settings <- data.frame(
  item = c(
    "score input", "clinical input", "feature set", "training fraction",
    "held-out test fraction", "matching", "training CV", "candidate model 1",
    "candidate model 2", "selection metric", "held-out metric", "DPS definition",
    "minimum group n", "seed", "input percentile scale"
  ),
  value = c(
    score_file,
    clinical_file,
    "24 bilateral hippocampal percentile features; age and sex excluded from predictors",
    training_fraction,
    test_fraction,
    "within each partition: exact sex plus nearest age, 1:1 HC, without replacement",
    paste0(cv_folds, "-fold cross-validation within the matched training set"),
    "regularized logistic regression (glmnet)",
    "single-hidden-layer neural network (nnet)",
    "highest training cross-validated ROC AUC",
    "ROC AUC from independent held-out probabilities",
    "held-out predicted probability of the target disease",
    minimum_group_n,
    random_seed,
    percentile_scale
  ),
  stringsAsFactors = FALSE
)
write.csv(settings, file.path(output_dir, "classification_settings.csv"), row.names = FALSE)

message(
  "Classification completed for ", nrow(metrics), " HC-versus-disease task(s). ",
  "Held-out DPS and metrics were written to: ", output_dir
)

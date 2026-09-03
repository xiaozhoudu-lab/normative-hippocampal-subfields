# ============================================================
# 04_Group_comparison.R
# Disease-versus-control comparison of regional normative
# percentiles and optional 24-ROI CMD values.
# ============================================================

rm(list = ls())

# Usage:
# Rscript 04_Group_comparison.R <normative_scores_long.csv> \
#   <clinical_data.csv> <output_dir> [individual_CMD.csv|-] \
#   [n_permutations] [seed] [diagnosis_column] [healthy_control_labels]
#
# The primary analysis uses 1:1 case-control matching with exact
# site/scanner and sex matching and nearest-age matching without
# replacement. Labels are permuted within matched pairs. Regional
# percentiles and CMD are kept in separate FDR families.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 04_Group_comparison.R <normative_scores_long.csv> ",
    "<clinical_data.csv> <output_dir> [individual_CMD.csv|-] ",
    "[n_permutations] [seed] [diagnosis_column] [healthy_control_labels]"
  )
}

score_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
cmd_arg <- if (length(args) >= 4) args[[4]] else "-"
n_permutations <- if (length(args) >= 5) as.integer(args[[5]]) else 5000L
random_seed <- if (length(args) >= 6) as.integer(args[[6]]) else 2025L
diagnosis_column <- if (length(args) >= 7) args[[7]] else "Diagnosis"
healthy_label_arg <- if (length(args) >= 8) args[[8]] else "HC"

if (!is.finite(n_permutations) || n_permutations < 1) {
  stop("n_permutations must be a positive integer")
}
if (!is.finite(random_seed)) stop("seed must be an integer")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
set.seed(random_seed)

minimum_total_n <- 10
minimum_group_n <- 2

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

make_site_id <- function(df) {
  existing <- if ("Site_ZZZ" %in% names(df)) safe_chr(df$Site_ZZZ) else rep("", nrow(df))
  if (all(c("Province", "Center", "Manufacturer") %in% names(df))) {
    components <- paste0(
      safe_chr(df$Province), safe_chr(df$Center), safe_chr(df$Manufacturer)
    )
  } else {
    components <- rep("", nrow(df))
  }
  ifelse(existing != "", existing, components)
}

normalize_diagnosis <- function(x) {
  out <- safe_chr(x)
  out[toupper(out) == "SVD"] <- "CSVD"
  out
}

format_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  if (length(x) == 1) return(sprintf("%.4f (NA)", x))
  sprintf("%.4f (%.4f)", mean(x), sd(x))
}

format_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  q <- quantile(x, c(0.25, 0.75), names = FALSE)
  sprintf("%.4f [%.4f-%.4f]", median(x), q[[1]], q[[2]])
}

welch_statistic <- function(y, group, orientation) {
  hc <- y[group == "HC"]
  disease <- y[group == "Disease"]
  hc <- hc[is.finite(hc)]
  disease <- disease[is.finite(disease)]
  if (length(hc) < 2 || length(disease) < 2) return(NA_real_)
  denominator <- sqrt(stats::var(hc) / length(hc) + stats::var(disease) / length(disease))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  difference <- if (orientation == "HC_minus_Disease") {
    mean(hc) - mean(disease)
  } else {
    mean(disease) - mean(hc)
  }
  difference / denominator
}

welch_details <- function(hc, disease, orientation) {
  hc <- hc[is.finite(hc)]
  disease <- disease[is.finite(disease)]
  v_hc <- stats::var(hc)
  v_disease <- stats::var(disease)
  se2 <- v_hc / length(hc) + v_disease / length(disease)
  difference <- if (orientation == "HC_minus_Disease") {
    mean(hc) - mean(disease)
  } else {
    mean(disease) - mean(hc)
  }
  statistic <- difference / sqrt(se2)
  degrees_freedom <- se2^2 / (
    (v_hc / length(hc))^2 / (length(hc) - 1) +
      (v_disease / length(disease))^2 / (length(disease) - 1)
  )
  p_value <- 2 * stats::pt(-abs(statistic), df = degrees_freedom)
  list(statistic = statistic, df = degrees_freedom, p_value = p_value)
}

cohens_d_details <- function(hc, disease, orientation) {
  hc <- hc[is.finite(hc)]
  disease <- disease[is.finite(disease)]
  n_hc <- length(hc)
  n_disease <- length(disease)
  pooled_variance <- (
    (n_hc - 1) * stats::var(hc) + (n_disease - 1) * stats::var(disease)
  ) / (n_hc + n_disease - 2)
  if (!is.finite(pooled_variance) || pooled_variance <= 0) {
    return(list(d = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  difference <- if (orientation == "HC_minus_Disease") {
    mean(hc) - mean(disease)
  } else {
    mean(disease) - mean(hc)
  }
  d <- difference / sqrt(pooled_variance)
  standard_error <- sqrt(
    (n_hc + n_disease) / (n_hc * n_disease) +
      d^2 / (2 * (n_hc + n_disease - 2))
  )
  list(d = d, lower = d - 1.96 * standard_error, upper = d + 1.96 * standard_error)
}

permute_within_pairs <- function(group, pair_id) {
  permuted <- group
  for (pair in unique(pair_id)) {
    index <- which(pair_id == pair)
    if (length(index) != 2L) stop("Every matched pair must contain exactly two participants")
    permuted[index] <- sample(group[index], replace = FALSE)
  }
  permuted
}

empty_retention <- function(disease, family, outcome, raw_hc, raw_disease, status) {
  data.frame(
    disease = disease,
    outcome_family = family,
    outcome = outcome,
    n_raw_HC = raw_hc,
    n_raw_Disease = raw_disease,
    n_matched_HC = 0L,
    n_matched_Disease = 0L,
    n_pairs = 0L,
    mean_absolute_age_difference = NA_real_,
    maximum_absolute_age_difference = NA_real_,
    status = status,
    stringsAsFactors = FALSE
  )
}

match_one_to_one <- function(data) {
  matched_rows <- list()
  pair_counter <- 0L
  exact_strata <- interaction(data$Site_ZZZ, data$Sex, drop = TRUE, lex.order = TRUE)

  for (stratum in levels(exact_strata)) {
    stratum_index <- which(exact_strata == stratum)
    case_index <- stratum_index[data$group[stratum_index] == "Disease"]
    control_index <- stratum_index[data$group[stratum_index] == "HC"]
    if (length(case_index) == 0L || length(control_index) == 0L) next

    case_index <- case_index[order(data$Age[case_index], data$match_id[case_index])]
    control_index <- control_index[order(data$Age[control_index], data$match_id[control_index])]

    while (length(case_index) > 0L && length(control_index) > 0L) {
      age_differences <- abs(outer(data$Age[case_index], data$Age[control_index], "-"))
      selected <- which(age_differences == min(age_differences), arr.ind = TRUE)[1, ]
      selected_case <- case_index[selected[[1]]]
      selected_control <- control_index[selected[[2]]]
      pair_counter <- pair_counter + 1L
      pair_id <- paste0("pair_", pair_counter)
      age_difference <- abs(data$Age[selected_case] - data$Age[selected_control])

      case_row <- data[selected_case, , drop = FALSE]
      control_row <- data[selected_control, , drop = FALSE]
      case_row$pair_id <- pair_id
      control_row$pair_id <- pair_id
      case_row$match_age_difference <- age_difference
      control_row$match_age_difference <- age_difference
      matched_rows[[length(matched_rows) + 1L]] <- case_row
      matched_rows[[length(matched_rows) + 1L]] <- control_row

      case_index <- case_index[-selected[[1]]]
      control_index <- control_index[-selected[[2]]]
    }
  }

  if (length(matched_rows) == 0L) return(data.frame())
  matched <- do.call(rbind, matched_rows)
  rownames(matched) <- NULL
  matched
}

analyze_matched <- function(data, disease, family, outcome, orientation) {
  sub <- data[
    data$Diagnosis %in% c(healthy_labels, disease) &
      is.finite(data$value) & is.finite(data$Age) &
      data$Sex != "" & data$Site_ZZZ != "",
    , drop = FALSE
  ]
  sub$group <- ifelse(sub$Diagnosis %in% healthy_labels, "HC", "Disease")
  raw_hc <- sum(sub$group == "HC")
  raw_disease <- sum(sub$group == "Disease")

  fail <- function(reason, status, retention = NULL) {
    if (is.null(retention)) {
      retention <- empty_retention(disease, family, outcome, raw_hc, raw_disease, status)
    }
    list(
      result = NULL,
      retention = retention,
      skip = data.frame(
        disease = disease, outcome_family = family, outcome = outcome,
        reason = reason, stringsAsFactors = FALSE
      )
    )
  }

  if (nrow(sub) < minimum_total_n || raw_hc < minimum_group_n || raw_disease < minimum_group_n) {
    return(fail("too few complete rows before 1:1 matching", "skipped_before_matching"))
  }

  sub_analysis <- match_one_to_one(sub)
  n_hc <- if (nrow(sub_analysis) > 0L) sum(sub_analysis$group == "HC") else 0L
  n_disease <- if (nrow(sub_analysis) > 0L) sum(sub_analysis$group == "Disease") else 0L
  n_pairs <- if (nrow(sub_analysis) > 0L) length(unique(sub_analysis$pair_id)) else 0L
  pair_age_differences <- if (nrow(sub_analysis) > 0L) {
    sub_analysis$match_age_difference[sub_analysis$group == "Disease"]
  } else {
    numeric(0)
  }

  retention <- data.frame(
    disease = disease,
    outcome_family = family,
    outcome = outcome,
    n_raw_HC = raw_hc,
    n_raw_Disease = raw_disease,
    n_matched_HC = n_hc,
    n_matched_Disease = n_disease,
    n_pairs = n_pairs,
    mean_absolute_age_difference = if (length(pair_age_differences)) {
      mean(pair_age_differences)
    } else {
      NA_real_
    },
    maximum_absolute_age_difference = if (length(pair_age_differences)) {
      max(pair_age_differences)
    } else {
      NA_real_
    },
    status = "pending",
    stringsAsFactors = FALSE
  )

  if (nrow(sub_analysis) < minimum_total_n ||
      n_hc < minimum_group_n || n_disease < minimum_group_n) {
    retention$status <- "skipped_after_matching"
    return(fail(
      "too few 1:1 pairs after exact site/sex and nearest-age matching",
      retention$status, retention
    ))
  }

  pair_sizes <- table(sub_analysis$pair_id)
  pair_composition <- table(sub_analysis$pair_id, sub_analysis$group)
  if (any(pair_sizes != 2L) ||
      any(pair_composition[, "HC"] != 1L) ||
      any(pair_composition[, "Disease"] != 1L)) {
    stop("Invalid 1:1 matched-pair construction for ", disease, " / ", outcome)
  }

  observed_t <- welch_statistic(sub_analysis$value, sub_analysis$group, orientation)
  if (!is.finite(observed_t)) {
    retention$status <- "skipped_nonfinite_statistic"
    return(fail("non-finite observed Welch statistic", retention$status, retention))
  }

  permutation_statistics <- numeric(n_permutations)
  for (iteration in seq_len(n_permutations)) {
    permuted_group <- permute_within_pairs(sub_analysis$group, sub_analysis$pair_id)
    permutation_statistics[[iteration]] <- welch_statistic(
      sub_analysis$value, permuted_group, orientation
    )
  }
  finite_permutations <- is.finite(permutation_statistics)
  empirical_p <- (
    sum(abs(permutation_statistics[finite_permutations]) >= abs(observed_t)) + 1
  ) / (sum(finite_permutations) + 1)

  hc <- sub_analysis$value[sub_analysis$group == "HC"]
  disease_values <- sub_analysis$value[sub_analysis$group == "Disease"]
  welch <- welch_details(hc, disease_values, orientation)
  effect <- cohens_d_details(hc, disease_values, orientation)
  calibration_status <- if ("calibration_status" %in% names(sub_analysis)) {
    safe_chr(sub_analysis$calibration_status)
  } else {
    rep("", nrow(sub_analysis))
  }
  uncalibrated <- grepl("uncalibrated", calibration_status, ignore.case = TRUE)

  retention$status <- "analyzed"
  result <- data.frame(
    disease = disease,
    outcome_family = family,
    outcome = outcome,
    comparison_direction = orientation,
    n_HC = n_hc,
    n_Disease = n_disease,
    n_pairs = n_pairs,
    mean_absolute_age_difference = mean(pair_age_differences),
    maximum_absolute_age_difference = max(pair_age_differences),
    n_uncalibrated_HC = sum(uncalibrated & sub_analysis$group == "HC"),
    n_uncalibrated_Disease = sum(uncalibrated & sub_analysis$group == "Disease"),
    mean_HC = mean(hc),
    mean_Disease = mean(disease_values),
    mean_sd_HC = format_mean_sd(hc),
    mean_sd_Disease = format_mean_sd(disease_values),
    median_iqr_HC = format_median_iqr(hc),
    median_iqr_Disease = format_median_iqr(disease_values),
    mean_difference_oriented = if (orientation == "HC_minus_Disease") {
      mean(hc) - mean(disease_values)
    } else {
      mean(disease_values) - mean(hc)
    },
    Welch_t = welch$statistic,
    Welch_df = welch$df,
    p_Welch = welch$p_value,
    p_permutation = empirical_p,
    Cohens_d = effect$d,
    Cohens_d_95CI_low = effect$lower,
    Cohens_d_95CI_high = effect$upper,
    stringsAsFactors = FALSE
  )
  list(result = result, retention = retention, skip = NULL)
}

clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
clinical$match_id <- make_match_id(clinical)
required_clinical_columns <- c("Age", "Sex", diagnosis_column)
missing_clinical_columns <- setdiff(required_clinical_columns, names(clinical))
if (length(missing_clinical_columns) > 0) {
  stop(
    "Clinical data are missing required columns: ",
    paste(missing_clinical_columns, collapse = ", ")
  )
}
clinical$Diagnosis <- normalize_diagnosis(clinical[[diagnosis_column]])
clinical$Age <- safe_num(clinical$Age)
clinical$Sex <- safe_chr(clinical$Sex)
clinical$Site_ZZZ <- make_site_id(clinical)
clinical <- clinical[clinical$match_id != "", , drop = FALSE]
if (anyDuplicated(clinical$match_id)) stop("Clinical data contain duplicate match_id values")
if (any(clinical$Site_ZZZ == "")) {
  stop("Some clinical rows have no Site_ZZZ or complete site components")
}
clinical <- clinical[, c("match_id", "Age", "Sex", "Site_ZZZ", "Diagnosis"), drop = FALSE]

healthy_labels <- trimws(unlist(strsplit(healthy_label_arg, ",", fixed = TRUE)))
healthy_labels <- normalize_diagnosis(healthy_labels[healthy_labels != ""])
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
scores$value <- safe_num(scores$Quant_score)
finite_percentiles <- scores$value[is.finite(scores$value)]
if (length(finite_percentiles) == 0) stop("Quant_score contains no finite values")
if (max(finite_percentiles) > 1.5) scores$value <- scores$value / 100
if (any(is.finite(scores$value) & (scores$value < 0 | scores$value > 1))) {
  stop("Quant_score contains values outside the detected percentile scale")
}
scores$outcome <- paste0(tolower(safe_chr(scores$hemisphere)), "_", safe_chr(scores$ROI))
duplicate_score <- duplicated(scores[c("match_id", "outcome")]) |
  duplicated(scores[c("match_id", "outcome")], fromLast = TRUE)
if (any(duplicate_score)) stop("Duplicate participant-region rows found in normative score input")

regional_columns <- c("match_id", "outcome", "value")
if ("calibration_status" %in% names(scores)) regional_columns <- c(regional_columns, "calibration_status")
regional_data <- merge(
  scores[, regional_columns, drop = FALSE], clinical,
  by = "match_id", all.x = TRUE, sort = FALSE
)
if (any(regional_data$Diagnosis == "" | is.na(regional_data$Diagnosis))) {
  warning("Some score rows could not be assigned a diagnosis and will not be analyzed")
}

diseases <- sort(unique(clinical$Diagnosis[
  clinical$Diagnosis != "" & !(clinical$Diagnosis %in% healthy_labels)
]))
if (length(diseases) == 0) stop("No disease labels remain after excluding healthy controls")

all_results <- list()
all_retention <- list()
all_skips <- list()

append_analysis <- function(analysis) {
  if (!is.null(analysis$result)) all_results[[length(all_results) + 1]] <<- analysis$result
  if (!is.null(analysis$retention)) all_retention[[length(all_retention) + 1]] <<- analysis$retention
  if (!is.null(analysis$skip)) all_skips[[length(all_skips) + 1]] <<- analysis$skip
}

regional_outcomes <- sort(unique(regional_data$outcome))
for (disease in diseases) {
  message("Regional comparison: ", disease)
  for (outcome in regional_outcomes) {
    outcome_data <- regional_data[regional_data$outcome == outcome, , drop = FALSE]
    append_analysis(analyze_matched(
      outcome_data, disease, "regional_percentile", outcome, "HC_minus_Disease"
    ))
  }
}

cmd_data <- NULL

if (!identical(cmd_arg, "-")) {
  cmd_file <- normalizePath(cmd_arg, winslash = "/", mustWork = TRUE)
  cmd_input <- read.csv(cmd_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("match_id", "CMD") %in% names(cmd_input))) {
    stop("CMD input must contain match_id and CMD columns from 05_CMD.R")
  }
  cmd_input$match_id <- safe_chr(cmd_input$match_id)
  if (anyDuplicated(cmd_input$match_id)) stop("CMD input contains duplicate match_id values")
  cmd_columns <- c("match_id", "CMD")
  if ("calibration_status" %in% names(cmd_input)) cmd_columns <- c(cmd_columns, "calibration_status")
  cmd_data <- merge(
    cmd_input[, cmd_columns, drop = FALSE], clinical,
    by = "match_id", all.x = TRUE, sort = FALSE
  )
  cmd_data$value <- safe_num(cmd_data$CMD)
  cmd_data$outcome <- "CMD_24ROI"

  for (disease in diseases) {
    message("CMD 1:1 matched comparison: ", disease)
    append_analysis(analyze_matched(
      cmd_data, disease, "CMD", "CMD_24ROI", "Disease_minus_HC"
    ))
  }
}

results <- if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
retention <- if (length(all_retention) > 0) do.call(rbind, all_retention) else data.frame()
skips <- if (length(all_skips) > 0) do.call(rbind, all_skips) else data.frame()

if (nrow(results) > 0) {
  results$pFDR <- NA_real_
  regional_index <- results$outcome_family == "regional_percentile"
  for (disease in unique(results$disease[regional_index])) {
    index <- which(regional_index & results$disease == disease)
    results$pFDR[index] <- stats::p.adjust(results$p_permutation[index], method = "BH")
  }
  cmd_index <- results$outcome_family == "CMD"
  if (any(cmd_index)) {
    results$pFDR[cmd_index] <- stats::p.adjust(
      results$p_permutation[cmd_index], method = "BH"
    )
  }
  results <- results[order(
    results$outcome_family, results$disease, results$pFDR, results$outcome
  ), , drop = FALSE]
  write.csv(results, file.path(output_dir, "group_comparison_all_outcomes.csv"), row.names = FALSE)
  write.csv(
    results[results$outcome_family == "regional_percentile", , drop = FALSE],
    file.path(output_dir, "regional_group_comparison.csv"), row.names = FALSE
  )
  if (any(results$outcome_family == "CMD")) {
    write.csv(
      results[results$outcome_family == "CMD", , drop = FALSE],
      file.path(output_dir, "CMD_group_comparison_1to1_matched.csv"), row.names = FALSE
    )
  }
}
if (nrow(retention) > 0) {
  write.csv(retention, file.path(output_dir, "group_comparison_sample_retention.csv"), row.names = FALSE)
}
if (nrow(skips) > 0) {
  write.csv(skips, file.path(output_dir, "group_comparison_skip_log.csv"), row.names = FALSE)
}
settings <- data.frame(
  item = c(
    "regional score input", "clinical input", "CMD input", "matching ratio",
    "exact matching", "age matching", "HC reuse", "age caliper",
    "permutation unit", "permutations", "seed",
    "regional effect direction", "CMD effect direction", "regional FDR family",
    "CMD FDR family"
  ),
  value = c(
    score_file, clinical_file, if (identical(cmd_arg, "-")) "not supplied" else cmd_arg,
    "1 case : 1 HC", "site/scanner and sex", "nearest absolute age",
    "without replacement", "none", "within matched pair",
    n_permutations, random_seed, "HC minus Disease", "Disease minus HC",
    "within each disease across regional outcomes",
    "across diseases for the CMD outcome"
  ),
  stringsAsFactors = FALSE
)
write.csv(settings, file.path(output_dir, "group_comparison_settings.csv"), row.names = FALSE)

message(
  "Group comparison completed. Regional results: ",
  sum(results$outcome_family == "regional_percentile"),
  "; 1:1 matched CMD results: ", sum(results$outcome_family == "CMD"),
  ". Outputs: ", output_dir
)

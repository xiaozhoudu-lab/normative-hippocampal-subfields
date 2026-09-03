# ============================================================
# 05_CMD.R
# Calculate 24-ROI centile Mahalanobis distance (CMD) from the
# percentile output of 02_Apply_normative_model.R.
# ============================================================

rm(list = ls())

# Usage:
# Rscript 05_CMD.R <normative_scores_long.csv> <clinical_data.csv|-> \
#   <output_dir> [diagnosis_column] [healthy_control_labels] [reference_rds|fit]
#
# Examples:
#   # Fit the CMD reference from complete healthy controls in this cohort.
#   Rscript 05_CMD.R scores/normative_scores_long.csv clinical.csv cmd_out
#
#   # Apply a previously saved CMD reference to new participants.
#   Rscript 05_CMD.R new_scores/normative_scores_long.csv new_clinical.csv \
#     new_cmd_out Diagnosis HC cmd_out/cmd_reference.rds
#
# Use "-" as clinical_data.csv when the score file already contains the
# diagnosis column, or when applying an existing reference without metadata.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 05_CMD.R <normative_scores_long.csv> ",
    "<clinical_data.csv|-> <output_dir> [diagnosis_column] ",
    "[healthy_control_labels] [reference_rds|fit]"
  )
}

score_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_arg <- args[[2]]
output_dir <- args[[3]]
diagnosis_column <- if (length(args) >= 4) args[[4]] else "Diagnosis"
healthy_labels <- if (length(args) >= 5) args[[5]] else "HC"
reference_arg <- if (length(args) >= 6) args[[6]] else "fit"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

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
  x <- sort(unique(safe_chr(x)))
  x <- x[x != ""]
  if (length(x) == 0) "" else paste(x, collapse = ";")
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

scores <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("match_id", "hemisphere", "ROI", "Quant_score")
missing_columns <- setdiff(required_columns, names(scores))
if (length(missing_columns) > 0) {
  stop(
    "The score file must be the long output from 02_Apply_normative_model.R; ",
    "missing columns: ", paste(missing_columns, collapse = ", ")
  )
}

scores$match_id <- safe_chr(scores$match_id)
scores$hemisphere_CMD <- normalize_hemisphere(scores$hemisphere)
scores$ROI_CMD <- normalize_roi(scores$ROI)
if (any(scores$match_id == "")) stop("The score file contains empty match_id values")
if (anyNA(scores$hemisphere_CMD)) {
  stop("Unrecognized hemisphere values: ", collapse_unique(scores$hemisphere[is.na(scores$hemisphere_CMD)]))
}
if (anyNA(scores$ROI_CMD)) {
  stop("Unrecognized ROI values: ", collapse_unique(scores$ROI[is.na(scores$ROI_CMD)]))
}
scores$feature_name_CMD <- paste0(scores$hemisphere_CMD, "_", scores$ROI_CMD)

unexpected_features <- setdiff(unique(scores$feature_name_CMD), feature_order)
if (length(unexpected_features) > 0) {
  stop("Unexpected CMD features: ", paste(unexpected_features, collapse = ", "))
}
duplicate_key <- duplicated(scores[c("match_id", "feature_name_CMD")]) |
  duplicated(scores[c("match_id", "feature_name_CMD")], fromLast = TRUE)
if (any(duplicate_key)) {
  duplicate_examples <- unique(paste(
    scores$match_id[duplicate_key], scores$feature_name_CMD[duplicate_key], sep = " / "
  ))
  stop(
    "Duplicate participant-feature rows found (values were not averaged): ",
    paste(head(duplicate_examples, 10), collapse = ", ")
  )
}

centile <- suppressWarnings(as.numeric(scores$Quant_score))
finite_centile <- centile[is.finite(centile)]
if (length(finite_centile) == 0) stop("Quant_score contains no finite values")
centile_scale <- if (max(finite_centile, na.rm = TRUE) > 1.5) "0_to_100" else "0_to_1"
if (centile_scale == "0_to_100") centile <- centile / 100
invalid_centile <- is.finite(centile) & (centile < 0 | centile > 1)
if (any(invalid_centile)) {
  stop("Quant_score contains values outside the detected percentile scale")
}
scores$CMD_centile <- centile

participant_ids <- unique(scores$match_id)
centile_matrix <- matrix(
  NA_real_, nrow = length(participant_ids), ncol = length(feature_order),
  dimnames = list(participant_ids, feature_order)
)
row_index <- match(scores$match_id, participant_ids)
column_index <- match(scores$feature_name_CMD, feature_order)
centile_matrix[cbind(row_index, column_index)] <- scores$CMD_centile

participant_data <- data.frame(
  match_id = participant_ids,
  n_CMD_features = rowSums(is.finite(centile_matrix)),
  complete_24ROI_profile = rowSums(is.finite(centile_matrix)) == length(feature_order),
  stringsAsFactors = FALSE
)
participant_data$missing_CMD_features <- vapply(
  seq_along(participant_ids),
  function(i) paste(feature_order[!is.finite(centile_matrix[i, ])], collapse = ";"),
  character(1)
)

if ("calibration_status" %in% names(scores)) {
  participant_data$calibration_status <- vapply(
    participant_ids,
    function(id) collapse_unique(scores$calibration_status[scores$match_id == id]),
    character(1)
  )
  participant_data$contains_uncalibrated_site_score <- grepl(
    "uncalibrated", participant_data$calibration_status, ignore.case = TRUE
  )
}

if (!identical(clinical_arg, "-")) {
  clinical_file <- normalizePath(clinical_arg, winslash = "/", mustWork = TRUE)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  clinical$match_id <- make_match_id(clinical)
  clinical <- clinical[clinical$match_id != "", , drop = FALSE]
  if (anyDuplicated(clinical$match_id)) {
    stop("Clinical data contain duplicate match_id values")
  }
  metadata_candidates <- unique(c(
    diagnosis_column, "Age", "Sex", "Site_ZZZ", "Province", "Center", "Manufacturer"
  ))
  metadata_columns <- intersect(metadata_candidates, names(clinical))
  clinical_keep <- clinical[, c("match_id", metadata_columns), drop = FALSE]
  participant_data$.input_order <- seq_len(nrow(participant_data))
  participant_data <- merge(
    participant_data, clinical_keep, by = "match_id", all.x = TRUE, sort = FALSE
  )
  participant_data <- participant_data[order(participant_data$.input_order), , drop = FALSE]
  participant_data$.input_order <- NULL
} else if (diagnosis_column %in% names(scores)) {
  diagnosis_by_id <- tapply(scores[[diagnosis_column]], scores$match_id, collapse_unique)
  participant_data[[diagnosis_column]] <- unname(diagnosis_by_id[participant_data$match_id])
}

complete_profile <- participant_data$complete_24ROI_profile
fit_reference <- tolower(reference_arg) %in% c("fit", "new", "none", "-")

if (fit_reference) {
  if (!(diagnosis_column %in% names(participant_data))) {
    stop(
      "Cannot fit a CMD reference without diagnosis labels. Supply a clinical CSV ",
      "containing column '", diagnosis_column, "', or load an existing reference RDS."
    )
  }
  hc_values <- trimws(unlist(strsplit(healthy_labels, ",", fixed = TRUE)))
  hc_values <- tolower(hc_values[hc_values != ""])
  if (length(hc_values) == 0) stop("At least one healthy-control label is required")
  is_healthy_control <- tolower(safe_chr(participant_data[[diagnosis_column]])) %in% hc_values
  reference_rows <- is_healthy_control & complete_profile
  n_reference <- sum(reference_rows)
  if (n_reference <= length(feature_order)) {
    stop(
      "Direct inversion of a 24-ROI covariance matrix requires more than 24 ",
      "complete healthy controls; found ", n_reference, "."
    )
  }

  healthy_matrix <- centile_matrix[
    participant_data$match_id[reference_rows], feature_order, drop = FALSE
  ]
  reference_mean <- colMeans(healthy_matrix)
  reference_covariance <- stats::cov(healthy_matrix)
  if (any(!is.finite(reference_covariance))) {
    stop("Healthy-control centile covariance contains non-finite values")
  }
  inversion_error <- NULL
  inverse_covariance <- tryCatch(
    solve(reference_covariance),
    error = function(error) {
      inversion_error <<- conditionMessage(error)
      NULL
    }
  )
  if (is.null(inverse_covariance) || any(!is.finite(inverse_covariance))) {
    stop(
      "The unregularized healthy-control centile covariance matrix could not be ",
      "inverted. This direct-centile CMD method cannot be calculated for the ",
      "current reference sample. Details: ",
      if (is.null(inversion_error)) "non-finite inverse" else inversion_error
    )
  }
  covariance_condition_number <- kappa(reference_covariance, exact = TRUE)
  cmd_reference <- list(
    version = "CMD_24ROI_direct_centile_v1",
    created = as.character(Sys.time()),
    feature_order = feature_order,
    transform = "none; direct 0-to-1 centile-score space",
    percentile_scale_in_fit_data = centile_scale,
    n_healthy_controls = n_reference,
    healthy_control_labels = hc_values,
    mean = reference_mean,
    covariance = reference_covariance,
    inverse_covariance = inverse_covariance,
    covariance_condition_number = covariance_condition_number
  )
  reference_path <- file.path(output_dir, "cmd_reference.rds")
  saveRDS(cmd_reference, reference_path)
  reference_source <- "fitted_from_current_healthy_controls"
} else {
  reference_path <- normalizePath(reference_arg, winslash = "/", mustWork = TRUE)
  cmd_reference <- readRDS(reference_path)
  required_reference_fields <- c(
    "version", "feature_order", "mean", "covariance", "inverse_covariance",
    "n_healthy_controls"
  )
  missing_reference_fields <- setdiff(required_reference_fields, names(cmd_reference))
  if (length(missing_reference_fields) > 0) {
    stop("CMD reference is missing fields: ", paste(missing_reference_fields, collapse = ", "))
  }
  if (!identical(as.character(cmd_reference$feature_order), feature_order)) {
    stop("CMD reference feature order does not match the required 24-ROI definition")
  }
  if (!identical(as.character(cmd_reference$version), "CMD_24ROI_direct_centile_v1")) {
    stop(
      "The supplied CMD reference was not fitted in direct centile-score space. ",
      "Refit the reference using the current 05_CMD.R script."
    )
  }
  reference_mean <- as.numeric(cmd_reference$mean)
  names(reference_mean) <- feature_order
  reference_covariance <- as.matrix(cmd_reference$covariance)
  inverse_covariance <- as.matrix(cmd_reference$inverse_covariance)
  if (length(reference_mean) != length(feature_order) ||
      !identical(dim(reference_covariance), c(length(feature_order), length(feature_order))) ||
      !identical(dim(inverse_covariance), c(length(feature_order), length(feature_order))) ||
      any(!is.finite(reference_mean)) || any(!is.finite(reference_covariance)) ||
      any(!is.finite(inverse_covariance))) {
    stop("CMD reference contains invalid mean or covariance dimensions/values")
  }
  covariance_condition_number <- if (!is.null(cmd_reference$covariance_condition_number)) {
    as.numeric(cmd_reference$covariance_condition_number)
  } else {
    kappa(reference_covariance, exact = TRUE)
  }
  n_reference <- as.integer(cmd_reference$n_healthy_controls)
  reference_source <- "loaded_external_reference"
  reference_rows <- rep(FALSE, nrow(participant_data))
}

cmd <- rep(NA_real_, nrow(participant_data))
complete_ids <- participant_data$match_id[complete_profile]
if (length(complete_ids) > 0) {
  x <- centile_matrix[complete_ids, feature_order, drop = FALSE]
  delta <- sweep(x, 2, reference_mean, FUN = "-")
  squared_distance <- rowSums((delta %*% inverse_covariance) * delta)
  squared_distance[squared_distance < 0 & squared_distance > -1e-8] <- 0
  if (any(squared_distance < 0, na.rm = TRUE)) {
    stop("Negative squared Mahalanobis distances indicate an invalid CMD reference")
  }
  cmd[match(complete_ids, participant_data$match_id)] <- sqrt(squared_distance)
}

participant_data$CMD <- cmd
participant_data$CMD_reference_member <- reference_rows
participant_data$CMD_reference_source <- reference_source
participant_data$CMD_reference_n_healthy_controls <- n_reference
participant_data$CMD_definition <- "24ROI_direct_centile_Mahalanobis"

write.csv(
  participant_data,
  file.path(output_dir, "individual_CMD.csv"),
  row.names = FALSE, na = ""
)

reference_mean_table <- data.frame(
  feature_name = feature_order,
  healthy_control_mean_centile = as.numeric(reference_mean),
  stringsAsFactors = FALSE
)
write.csv(
  reference_mean_table,
  file.path(output_dir, "cmd_reference_mean.csv"),
  row.names = FALSE
)

covariance_table <- data.frame(
  feature_name = feature_order,
  as.data.frame(reference_covariance, check.names = FALSE),
  check.names = FALSE
)
names(covariance_table)[-1] <- feature_order
write.csv(
  covariance_table,
  file.path(output_dir, "cmd_reference_covariance.csv"),
  row.names = FALSE
)

qa <- data.frame(
  item = c(
    "CMD definition", "ROI vector", "analysis space", "input percentile scale",
    "complete profiles", "incomplete profiles", "reference source",
    "healthy controls in reference", "covariance estimator", "matrix inverse",
    "covariance condition number", "reference file"
  ),
  value = c(
    "sqrt((x - mean_HC)' inverse(cov_HC) (x - mean_HC))",
    paste(length(feature_order), "features: bilateral whole hippocampus plus 11 subfields"),
    "direct 0-to-1 centile scores; no transformation or clipping",
    centile_scale,
    sum(complete_profile),
    sum(!complete_profile),
    reference_source,
    n_reference,
    "unregularized sample covariance among complete healthy controls",
    "direct solve(covariance); no ridge and no pseudoinverse",
    covariance_condition_number,
    reference_path
  ),
  stringsAsFactors = FALSE
)
write.csv(qa, file.path(output_dir, "cmd_QA.csv"), row.names = FALSE)

message(
  "CMD completed for ", sum(is.finite(participant_data$CMD)), " of ",
  nrow(participant_data), " participant(s). Outputs: ", output_dir
)

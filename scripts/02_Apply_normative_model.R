# ============================================================
# 02_Apply_normative_model.R
# Apply fitted hippocampal normative models to one participant
# or a cohort and calculate auditable Z and percentile scores.
# ============================================================

rm(list = ls())

# Usage:
# Rscript 02_Apply_normative_model.R <source_dir> <clinical_data.csv> \
#   <mr_workbook.xlsx> <model_root> <output_dir> [primary|icv]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 02_Apply_normative_model.R <source_dir> ",
    "<clinical_data.csv> <mr_workbook.xlsx> <model_root> <output_dir> ",
    "[primary|icv]"
  )
}

source_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
mr_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
model_root <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
output_dir <- args[[5]]
model_variant <- if (length(args) >= 6) tolower(args[[6]]) else "primary"
if (!(model_variant %in% c("primary", "icv"))) {
  stop("model variant must be 'primary' or 'icv'")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(gamlss)
  library(gamlss.add)
})

setwd(source_dir)
source_files <- c(
  "100.common-variables.r", "101.common-functions.r", "102.gamlss-recode.r",
  "300.variables.r", "301.functions.r", "ZZZ_function.R"
)
for (ff in source_files) {
  if (!file.exists(ff)) stop("Required helper file not found: ", file.path(source_dir, ff))
  source(ff)
}

sheets <- c("lh.hipposubfields.vol.table", "rh.hipposubfields.vol.table")
target_features <- c(
  "subiculum", "CA1", "presubiculum", "molecular_layer_HP",
  "GC.ML.DG", "CA3", "CA4", "Hippocampal_tail",
  "parasubiculum", "fimbria", "HATA", "Whole_hippocampus"
)

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

make_match_id <- function(df, clinical = FALSE) {
  candidates <- if (clinical) {
    list(
      c("Freesufer_Path2", "Freesufer_Path3"),
      c("Freesurfer_Path2", "Freesurfer_Path3")
    )
  } else {
    list(c("Freesurfer_Path2", "Freesurfer_Path3"))
  }
  for (pair in candidates) {
    if (all(pair %in% names(df))) {
      return(paste0(safe_chr(df[[pair[1]]]), safe_chr(df[[pair[2]]])))
    }
  }
  stop("Cannot identify FreeSurfer path columns used to match clinical and MRI data")
}

make_site_id <- function(df) {
  existing <- if ("Site_ZZZ" %in% names(df)) safe_chr(df$Site_ZZZ) else rep("", nrow(df))
  if (all(c("Province", "Center", "Manufacturer") %in% names(df))) {
    from_components <- paste0(
      safe_chr(df$Province), safe_chr(df$Center), safe_chr(df$Manufacturer)
    )
  } else {
    from_components <- rep("", nrow(df))
  }
  site <- ifelse(existing != "", existing, from_components)
  if (any(site == "")) {
    stop("Some rows have neither Site_ZZZ nor complete Province/Center/Manufacturer fields")
  }
  site
}

get_site_coef <- function(model1, what_name, site_levels) {
  smo <- model1[[paste0(what_name, ".coefSmo")]]
  if (is.null(smo) || length(smo) == 0) return(NULL)
  for (k in seq_along(smo)) {
    cc <- smo[[k]]$coef
    if (is.null(cc) || is.null(names(cc))) next
    hit <- intersect(names(cc), site_levels)
    if (length(hit) > 0) return(cc[hit])
  }
  NULL
}

choose_neutral_site <- function(model1, training_data) {
  site_levels <- levels(training_data$Site_ZZZ)
  if (is.null(site_levels) || length(site_levels) == 0) {
    stop("Training data contain no Site_ZZZ factor levels")
  }
  mu_coef <- get_site_coef(model1, "mu", site_levels)
  sigma_coef <- get_site_coef(model1, "sigma", site_levels)
  neutral_site <- NA_character_
  method <- "fallback_most_common_training_site"

  if (!is.null(mu_coef) && !is.null(sigma_coef)) {
    common_sites <- intersect(names(mu_coef), names(sigma_coef))
    if (length(common_sites) > 0) {
      neutral_score <-
        abs(mu_coef[common_sites] - mean(mu_coef[common_sites], na.rm = TRUE)) +
        abs(sigma_coef[common_sites] - mean(sigma_coef[common_sites], na.rm = TRUE))
      neutral_site <- common_sites[which.min(neutral_score)]
      method <- "mu_sigma_closest_to_mean_random_effect"
    }
  }
  if (is.na(neutral_site) && !is.null(mu_coef) && length(mu_coef) > 0) {
    neutral_site <- names(mu_coef)[which.min(abs(mu_coef - mean(mu_coef, na.rm = TRUE)))]
    method <- "mu_closest_to_mean_random_effect"
  }
  if (is.na(neutral_site) && !is.null(sigma_coef) && length(sigma_coef) > 0) {
    neutral_site <- names(sigma_coef)[which.min(abs(sigma_coef - mean(sigma_coef, na.rm = TRUE)))]
    method <- "sigma_closest_to_mean_random_effect"
  }
  if (is.na(neutral_site)) {
    site_tab <- sort(table(as.character(training_data$Site_ZZZ)), decreasing = TRUE)
    neutral_site <- if (length(site_tab) > 0) names(site_tab)[1] else site_levels[1]
  }
  list(neutral_site = neutral_site, method = method)
}

reset_gamlss_formula_environments <- function(model1, env = .GlobalEnv) {
  formula_slots <- names(model1)[vapply(model1, inherits, logical(1), what = "formula")]
  for (slot in formula_slots) environment(model1[[slot]]) <- env
  term_slots <- names(model1)[vapply(model1, inherits, logical(1), what = "terms")]
  for (slot in term_slots) attr(model1[[slot]], ".Environment") <- env
  model1
}

prime_gamlss_prediction_environment <- function(model1, training_data, m0 = NULL) {
  assign("data1", training_data, envir = .GlobalEnv)
  assign("object", model1, envir = .GlobalEnv)
  if (!is.null(m0)) assign("m0", m0, envir = .GlobalEnv)
  for (pkg_env in list(asNamespace("gamlss"), .GlobalEnv)) {
    for (env_name in c("gamlss.environment", "gamlss.env")) {
      if (exists(env_name, envir = pkg_env, inherits = FALSE)) {
        pred_env <- get(env_name, envir = pkg_env, inherits = FALSE)
        if (is.environment(pred_env)) {
          try(assign("object", model1, envir = pred_env), silent = TRUE)
          try(assign("data1", training_data, envir = pred_env), silent = TRUE)
          if (!is.null(m0)) try(assign("m0", m0, envir = pred_env), silent = TRUE)
        }
      }
    }
  }
  invisible(TRUE)
}

predict_gamlss_response <- function(model1, newdata, training_data, what, m0 = NULL) {
  model1 <- reset_gamlss_formula_environments(model1)
  prime_gamlss_prediction_environment(model1, training_data, m0)
  extra_cols <- setdiff(names(newdata), names(training_data))
  if (length(extra_cols) > 0) {
    stop("Prediction data contain columns absent from training data: ", paste(extra_cols, collapse = ", "))
  }
  as.numeric(predict(model1, newdata = newdata, type = "response", what = what))
}

make_model_newdata <- function(df, training_data, site_used) {
  required <- c("Age", "Sex", "Site_ZZZ", "tem_feature")
  if ("ICV" %in% names(training_data)) required <- c(required, "ICV")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Scoring data missing model covariates: ", paste(missing_cols, collapse = ", "))
  }
  nd <- df[, required, drop = FALSE]
  nd$Site_ZZZ <- factor(site_used, levels = levels(training_data$Site_ZZZ))
  nd$Sex <- factor(as.character(nd$Sex), levels = levels(training_data$Sex))
  if (anyNA(nd$Site_ZZZ)) stop("Prediction site is not available in model factor levels")
  if (anyNA(nd$Sex)) stop("Scoring data contain Sex values absent from model training data")
  nd$feature <- nd$tem_feature
  keep <- names(nd)[names(nd) %in% names(training_data)]
  nd[, keep, drop = FALSE]
}

calculate_scores <- function(model1, y, age, mu, sigma, nu) {
  z_score <- zzz_cent(
    obj = model1, type = "z-scores", mu = mu, sigma = sigma, nu = nu,
    xname = "Age", xvalues = age, yval = y,
    calibration = FALSE, lpar = 3
  )
  percentile <- zzz_cent(
    obj = model1, type = "z-scores", mu = mu, sigma = sigma, nu = nu,
    xname = "Age", xvalues = age, yval = y,
    calibration = FALSE, lpar = 3, cdf = TRUE
  )
  data.frame(Z_score = as.numeric(z_score), Quant_score = as.numeric(percentile))
}

find_model_file <- function(sheet, roi) {
  calibrated <- file.path(
    model_root, sheet,
    paste0(sheet, "_", roi, "_new_site_calibrated_results.rds")
  )
  base <- if (model_variant == "icv") {
    file.path(model_root, sheet, paste0(sheet, "_", roi, "_loop_our_model_ICV.rds"))
  } else {
    file.path(model_root, sheet, paste0(sheet, "_", roi, "_loop_our_model.rds"))
  }
  if (file.exists(calibrated)) return(calibrated)
  if (file.exists(base)) return(base)
  stop("Model file not found for ", sheet, " / ", roi, " under ", model_root)
}

clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
clinical$match_id <- make_match_id(clinical, clinical = TRUE)
clinical$Site_ZZZ <- make_site_id(clinical)
clinical$Age <- safe_num(clinical$Age)
clinical$Sex <- safe_chr(clinical$Sex)
if (!("ICV" %in% names(clinical)) && "EstimatedTotalIntraCranialVol" %in% names(clinical)) {
  clinical$ICV <- safe_num(clinical$EstimatedTotalIntraCranialVol)
}
clinical <- clinical[clinical$match_id != "", , drop = FALSE]
clinical <- clinical[!duplicated(clinical$match_id), , drop = FALSE]
rownames(clinical) <- clinical$match_id

all_scores <- list()

for (sheet in sheets) {
  message("Processing sheet: ", sheet)
  MRI <- as.data.frame(read_excel(mr_file, sheet = sheet))
  MRI$match_id <- make_match_id(MRI, clinical = FALSE)
  MRI <- MRI[MRI$match_id != "", , drop = FALSE]
  MRI <- MRI[!duplicated(MRI$match_id), , drop = FALSE]
  rownames(MRI) <- MRI$match_id

  missing_features <- setdiff(target_features, names(MRI))
  if (length(missing_features) > 0) {
    stop("Missing hippocampal features in ", sheet, ": ", paste(missing_features, collapse = ", "))
  }
  matched_ids <- intersect(rownames(clinical), rownames(MRI))
  if (length(matched_ids) == 0) stop("No matched clinical/MRI rows for ", sheet)

  for (roi in target_features) {
    model_file <- find_model_file(sheet, roi)
    results <- readRDS(model_file)
    if (is.null(results$m2) || is.null(results$data1)) {
      stop("Model RDS lacks results$m2 or results$data1: ", model_file)
    }
    model_is_icv <- isTRUE(results$ICV_adjusted) || "ICV" %in% names(results$data1)
    if ((model_variant == "icv") != model_is_icv) {
      stop("Model variant mismatch for ", model_file)
    }

    model1 <- results$m2
    m0 <- results$m0
    training_data <- as.data.frame(results$data1)
    training_data$Site_ZZZ <- droplevels(factor(as.character(training_data$Site_ZZZ)))
    training_data$Sex <- droplevels(factor(as.character(training_data$Sex)))

    scoring_data <- clinical[matched_ids, , drop = FALSE]
    scoring_data$tem_feature <- safe_num(MRI[matched_ids, roi])
    valid <- is.finite(scoring_data$Age) & is.finite(scoring_data$tem_feature) &
      scoring_data$Sex != "" & scoring_data$Site_ZZZ != ""
    if (model_is_icv) {
      if (!("ICV" %in% names(scoring_data))) {
        stop("ICV model requested but scoring data contain no ICV column")
      }
      scoring_data$ICV <- safe_num(scoring_data$ICV)
      valid <- valid & is.finite(scoring_data$ICV) & scoring_data$ICV > 0
    }
    scoring_data <- scoring_data[valid, , drop = FALSE]
    if (nrow(scoring_data) == 0) next

    neutral <- choose_neutral_site(model1, training_data)
    site_original <- safe_chr(scoring_data$Site_ZZZ)
    site_seen <- site_original %in% levels(training_data$Site_ZZZ)
    site_used <- ifelse(site_seen, site_original, neutral$neutral_site)

    newdata <- make_model_newdata(scoring_data, training_data, site_used)
    mu <- predict_gamlss_response(model1, newdata, training_data, "mu", m0)
    sigma <- predict_gamlss_response(model1, newdata, training_data, "sigma", m0)
    nu <- predict_gamlss_response(model1, newdata, training_data, "nu", m0)

    mu_offset <- rep(0, nrow(scoring_data))
    sigma_offset <- rep(0, nrow(scoring_data))
    calibration_n_hc <- integer(nrow(scoring_data))
    status <- ifelse(site_seen, "seen_training_site", "uncalibrated_neutral_fallback")

    calibration_table <- results$site_calibration
    if (!is.null(calibration_table) && nrow(calibration_table) > 0) {
      for (new_site in unique(site_original[!site_seen])) {
        hit <- calibration_table[
          calibration_table$Site_ZZZ == new_site &
            calibration_table$status == "calibrated_from_healthy_controls",
          , drop = FALSE
        ]
        if (nrow(hit) == 0) next
        if (nrow(hit) > 1) stop("Duplicate calibration rows for site ", new_site, " / ", roi)
        idx <- which(site_original == new_site)
        mu_hit <- safe_num(hit$mu_log_offset[1])
        sigma_hit <- safe_num(hit$sigma_log_offset[1])
        if (!is.finite(mu_hit) || !is.finite(sigma_hit)) {
          stop("Non-finite calibration offsets for site ", new_site, " / ", roi)
        }
        mu_offset[idx] <- mu_hit
        sigma_offset[idx] <- sigma_hit
        calibration_n_hc[idx] <- as.integer(hit$n_healthy_controls[1])
        status[idx] <- "calibrated_from_healthy_controls"
      }
    }

    mu <- mu * exp(mu_offset)
    sigma <- sigma * exp(sigma_offset)
    scores <- calculate_scores(
      model1, scoring_data$tem_feature, scoring_data$Age, mu, sigma, nu
    )

    hemisphere <- if (startsWith(sheet, "lh.")) "left" else "right"
    score_part <- data.frame(
      match_id = rownames(scoring_data),
      hemisphere = hemisphere,
      ROI = roi,
      Age = scoring_data$Age,
      Sex = scoring_data$Sex,
      Site_ZZZ_original = site_original,
      Site_ZZZ_used_for_base_prediction = site_used,
      Site_seen_in_training = site_seen,
      calibration_status = status,
      calibration_n_healthy_controls = calibration_n_hc,
      neutral_site_used = neutral$neutral_site,
      neutral_site_method = neutral$method,
      volume = scoring_data$tem_feature,
      Z_score = scores$Z_score,
      Quant_score = scores$Quant_score,
      stringsAsFactors = FALSE
    )
    score_part$feature_name <- paste0(hemisphere, "_", roi)
    all_scores[[paste(sheet, roi, sep = "__")]] <- score_part
  }
}

scores_long <- bind_rows(all_scores)
if (nrow(scores_long) == 0) stop("No participants were scored")
write.csv(scores_long, file.path(output_dir, "normative_scores_long.csv"), row.names = FALSE)

percentiles_wide <- reshape(
  scores_long[, c("match_id", "feature_name", "Quant_score")],
  idvar = "match_id", timevar = "feature_name", direction = "wide"
)
names(percentiles_wide) <- sub("^Quant_score\\.", "", names(percentiles_wide))
write.csv(percentiles_wide, file.path(output_dir, "normative_percentiles_wide.csv"), row.names = FALSE)

zscores_wide <- reshape(
  scores_long[, c("match_id", "feature_name", "Z_score")],
  idvar = "match_id", timevar = "feature_name", direction = "wide"
)
names(zscores_wide) <- sub("^Z_score\\.", "", names(zscores_wide))
write.csv(zscores_wide, file.path(output_dir, "normative_zscores_wide.csv"), row.names = FALSE)

audit <- scores_long[, c(
  "match_id", "hemisphere", "ROI", "Site_ZZZ_original",
  "Site_ZZZ_used_for_base_prediction", "Site_seen_in_training",
  "calibration_status", "calibration_n_healthy_controls",
  "neutral_site_used", "neutral_site_method"
)]
write.csv(audit, file.path(output_dir, "site_scoring_audit.csv"), row.names = FALSE)

message("Normative scoring completed for ", length(unique(scores_long$match_id)), " participant(s).")

# ============================================================
# 03_New_site_calibration.R
# Calibrate hippocampal normative models for unseen sites and
# calculate auditable Z scores and percentile scores.
# ============================================================

rm(list = ls())

# Usage:
# Rscript 03_New_site_calibration.R <source_dir> <clinical_data.csv> \
#   <mr_workbook.xlsx> <model_root> <output_dir> \
#   [diagnosis_column] [healthy_labels_csv] [minimum_controls]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 03_New_site_calibration.R <source_dir> ",
    "<clinical_data.csv> <mr_workbook.xlsx> <model_root> <output_dir> ",
    "[diagnosis_column] [healthy_labels_csv] [minimum_controls]"
  )
}

source_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
mr_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
model_root <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
output_root <- args[[5]]
diagnosis_column <- if (length(args) >= 6) args[[6]] else "Diagnosis"
healthy_labels <- if (length(args) >= 7) {
  trimws(strsplit(args[[7]], ",", fixed = TRUE)[[1]])
} else {
  c("HC", "NC", "CN", "Control", "Healthy Control", "Healthy_Control")
}
minimum_controls <- if (length(args) >= 8) as.integer(args[[8]]) else 20L
if (is.na(minimum_controls) || minimum_controls < 5) {
  stop("minimum_controls must be an integer of at least 5")
}
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(gamlss)
  library(gamlss.add)
  library(stringr)
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

normalise_label <- function(x) {
  tolower(gsub("[[:space:]_-]+", "", safe_chr(x)))
}

make_site_id <- function(df) {
  existing <- if ("Site_ZZZ" %in% names(df)) safe_chr(df$Site_ZZZ) else rep("", nrow(df))
  component_names <- c("Province", "Center", "Manufacturer")
  if (all(component_names %in% names(df))) {
    from_components <- paste0(
      safe_chr(df$Province), safe_chr(df$Center), safe_chr(df$Manufacturer)
    )
  } else {
    from_components <- rep("", nrow(df))
  }
  site <- ifelse(existing != "", existing, from_components)
  if (any(site == "")) {
    stop("Some rows have neither Site_ZZZ nor complete Province/Center/Manufacturer site information")
  }
  site
}

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
  stop("Cannot identify FreeSurfer path columns needed to match clinical and MRI data")
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
  if (anyNA(nd$Site_ZZZ)) stop("Internal error: prediction site is not a training factor level")
  if (anyNA(nd$Sex)) stop("Scoring data contain Sex values absent from model training data")
  nd$feature <- nd$tem_feature
  keep <- names(nd)[names(nd) %in% names(training_data)]
  nd[, keep, drop = FALSE]
}

estimate_site_offsets <- function(y, mu_base, sigma_base, nu_base) {
  ok <- is.finite(y) & is.finite(mu_base) & is.finite(sigma_base) &
    is.finite(nu_base) & mu_base > 0 & sigma_base > 0
  y <- y[ok]
  mu_base <- mu_base[ok]
  sigma_base <- sigma_base[ok]
  nu_base <- nu_base[ok]
  if (length(y) < 5) return(list(success = FALSE, message = "too_few_valid_controls"))

  objective <- function(delta) {
    density <- suppressWarnings(gamlss.dist::dGG(
      y,
      mu = mu_base * exp(delta[1]),
      sigma = sigma_base * exp(delta[2]),
      nu = nu_base,
      log = TRUE
    ))
    if (any(!is.finite(density))) return(.Machine$double.xmax / 100)
    -sum(density) + 0.01 * sum(delta^2)
  }

  fit <- try(
    optim(c(mu = 0, sigma = 0), objective, method = "BFGS", control = list(maxit = 500)),
    silent = TRUE
  )
  if (inherits(fit, "try-error") || fit$convergence != 0 || any(!is.finite(fit$par))) {
    return(list(success = FALSE, message = "optimisation_failed"))
  }
  list(
    success = TRUE,
    mu_log_offset = unname(fit$par[1]),
    sigma_log_offset = unname(fit$par[2]),
    objective = fit$value,
    message = "calibrated_from_healthy_controls"
  )
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
  data.frame(
    Z_score = as.numeric(z_score),
    Quant_score = as.numeric(percentile)
  )
}

clinical_all <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!(diagnosis_column %in% names(clinical_all))) {
  stop(
    "Diagnosis column not found: ", diagnosis_column,
    ". It is required so only healthy controls are used for calibration."
  )
}
clinical_all$match_id <- make_match_id(clinical_all, clinical = TRUE)
clinical_all$Site_ZZZ <- make_site_id(clinical_all)
clinical_all$Age <- safe_num(clinical_all$Age)
clinical_all$Sex <- safe_chr(clinical_all$Sex)
clinical_all$.Diagnosis_for_calibration <- safe_chr(clinical_all[[diagnosis_column]])
clinical_all$.is_healthy_control <-
  normalise_label(clinical_all$.Diagnosis_for_calibration) %in% normalise_label(healthy_labels)
clinical_all <- clinical_all[clinical_all$match_id != "", , drop = FALSE]
clinical_all <- clinical_all[!duplicated(clinical_all$match_id), , drop = FALSE]
rownames(clinical_all) <- clinical_all$match_id

all_calibration_tables <- list()
all_qa_tables <- list()

for (sheet in sheets) {
  message("Processing sheet: ", sheet)
  sheet_output <- file.path(output_root, sheet)
  if (!dir.exists(sheet_output)) dir.create(sheet_output, recursive = TRUE)

  MRI <- as.data.frame(read_excel(mr_file, sheet = sheet))
  MRI$match_id <- make_match_id(MRI, clinical = FALSE)
  MRI <- MRI[MRI$match_id != "", , drop = FALSE]
  MRI <- MRI[!duplicated(MRI$match_id), , drop = FALSE]
  rownames(MRI) <- MRI$match_id
  missing_features <- setdiff(target_features, names(MRI))
  if (length(missing_features) > 0) {
    stop("Missing hippocampal features in ", sheet, ": ", paste(missing_features, collapse = ", "))
  }
  matched_ids <- intersect(rownames(clinical_all), rownames(MRI))
  if (length(matched_ids) == 0) stop("No matched clinical/MRI rows for ", sheet)

  for (roi in target_features) {
    message("  ROI: ", roi)
    model_candidates <- c(
      file.path(model_root, sheet, paste0(sheet, "_", roi, "_loop_our_model.rds")),
      file.path(model_root, sheet, paste0(sheet, "_", roi, "_loop_our_model_ICV.rds"))
    )
    existing_models <- model_candidates[file.exists(model_candidates)]
    if (length(existing_models) == 0) stop("Model file not found for ", sheet, " / ", roi)
    if (length(existing_models) > 1) {
      stop(
        "Both primary and ICV-adjusted models were found for ", sheet, " / ", roi,
        ". Use a model_root containing only the intended model family."
      )
    }
    model_file <- existing_models[1]

    results <- readRDS(model_file)
    if (is.null(results$m2) || is.null(results$data1)) {
      stop("Model RDS lacks results$m2 or results$data1: ", model_file)
    }
    model1 <- results$m2
    m0 <- results$m0
    training_data <- as.data.frame(results$data1)
    training_data$Site_ZZZ <- droplevels(factor(as.character(training_data$Site_ZZZ)))
    training_data$Sex <- droplevels(factor(as.character(training_data$Sex)))

    scoring_data <- clinical_all[matched_ids, , drop = FALSE]
    scoring_data$tem_feature <- safe_num(MRI[matched_ids, roi])
    scoring_data <- scoring_data[
      is.finite(scoring_data$Age) & is.finite(scoring_data$tem_feature) &
        scoring_data$Sex != "" & scoring_data$Site_ZZZ != "",
      , drop = FALSE
    ]
    if ("ICV" %in% names(training_data)) {
      if (!("ICV" %in% names(scoring_data))) {
        stop("ICV-adjusted model selected but clinical scoring data contain no ICV column")
      }
      scoring_data$ICV <- safe_num(scoring_data$ICV)
      scoring_data <- scoring_data[is.finite(scoring_data$ICV) & scoring_data$ICV > 0, , drop = FALSE]
    }

    neutral <- choose_neutral_site(model1, training_data)
    site_original <- safe_chr(scoring_data$Site_ZZZ)
    site_seen <- site_original %in% levels(training_data$Site_ZZZ)
    site_used_for_base <- ifelse(site_seen, site_original, neutral$neutral_site)
    newdata <- make_model_newdata(scoring_data, training_data, site_used_for_base)
    mu_base <- predict_gamlss_response(model1, newdata, training_data, "mu", m0)
    sigma_base <- predict_gamlss_response(model1, newdata, training_data, "sigma", m0)
    nu_base <- predict_gamlss_response(model1, newdata, training_data, "nu", m0)

    mu_offset <- rep(0, nrow(scoring_data))
    sigma_offset <- rep(0, nrow(scoring_data))
    calibration_status <- ifelse(site_seen, "seen_training_site", "uncalibrated_neutral_fallback")
    calibration_n_hc <- integer(nrow(scoring_data))
    calibration_rows <- list()

    new_sites <- sort(unique(site_original[!site_seen]))
    for (new_site in new_sites) {
      idx_site <- which(site_original == new_site)
      idx_hc <- idx_site[scoring_data$.is_healthy_control[idx_site]]
      n_hc <- length(idx_hc)
      calibration_n_hc[idx_site] <- n_hc

      if (n_hc < minimum_controls) {
        warning(
          "Site ", new_site, " has only ", n_hc,
          " healthy controls; using explicit neutral-site fallback for ", sheet, " / ", roi
        )
        calibration_rows[[new_site]] <- data.frame(
          Site_ZZZ = new_site, n_healthy_controls = n_hc,
          mu_log_offset = 0, sigma_log_offset = 0,
          status = "insufficient_controls_neutral_fallback",
          neutral_site = neutral$neutral_site, neutral_site_method = neutral$method,
          stringsAsFactors = FALSE
        )
        next
      }

      fit <- estimate_site_offsets(
        scoring_data$tem_feature[idx_hc],
        mu_base[idx_hc], sigma_base[idx_hc], nu_base[idx_hc]
      )
      if (!isTRUE(fit$success)) {
        warning("Calibration failed for site ", new_site, "; using neutral-site fallback")
        calibration_rows[[new_site]] <- data.frame(
          Site_ZZZ = new_site, n_healthy_controls = n_hc,
          mu_log_offset = 0, sigma_log_offset = 0, status = fit$message,
          neutral_site = neutral$neutral_site, neutral_site_method = neutral$method,
          stringsAsFactors = FALSE
        )
        next
      }

      mu_offset[idx_site] <- fit$mu_log_offset
      sigma_offset[idx_site] <- fit$sigma_log_offset
      calibration_status[idx_site] <- "calibrated_from_healthy_controls"
      calibration_rows[[new_site]] <- data.frame(
        Site_ZZZ = new_site, n_healthy_controls = n_hc,
        mu_log_offset = fit$mu_log_offset, sigma_log_offset = fit$sigma_log_offset,
        status = fit$message, neutral_site = neutral$neutral_site,
        neutral_site_method = neutral$method, stringsAsFactors = FALSE
      )
    }

    mu_calibrated <- mu_base * exp(mu_offset)
    sigma_calibrated <- sigma_base * exp(sigma_offset)
    scores <- calculate_scores(
      model1,
      scoring_data$tem_feature,
      scoring_data$Age,
      mu_calibrated, sigma_calibrated, nu_base
    )

    score_output <- scoring_data
    score_output$Site_ZZZ_original <- site_original
    score_output$Site_ZZZ_used_for_base_prediction <- site_used_for_base
    score_output$Site_seen_in_training <- site_seen
    score_output$calibration_status <- calibration_status
    score_output$calibration_n_healthy_controls <- calibration_n_hc
    score_output$mu_log_offset <- mu_offset
    score_output$sigma_log_offset <- sigma_offset
    score_output$mu_calibrated <- mu_calibrated
    score_output$sigma_calibrated <- sigma_calibrated
    score_output$nu_calibrated <- nu_base
    score_output$Z_score <- scores$Z_score
    score_output$Quant_score <- scores$Quant_score

    calibration_table <- if (length(calibration_rows) > 0) {
      bind_rows(calibration_rows)
    } else {
      data.frame(
        Site_ZZZ = character(), n_healthy_controls = integer(),
        mu_log_offset = numeric(), sigma_log_offset = numeric(), status = character(),
        neutral_site = character(), neutral_site_method = character(),
        stringsAsFactors = FALSE
      )
    }
    calibration_table$Sheet <- sheet
    calibration_table$ROI <- roi

    qa_table <- score_output %>%
      filter(.is_healthy_control) %>%
      group_by(Site_ZZZ_original, calibration_status) %>%
      summarise(
        n_healthy_controls = n(),
        median_percentile = median(Quant_score, na.rm = TRUE),
        q1_percentile = quantile(Quant_score, 0.25, na.rm = TRUE),
        q3_percentile = quantile(Quant_score, 0.75, na.rm = TRUE),
        pct_below_0.025 = mean(Quant_score < 0.025, na.rm = TRUE),
        pct_above_0.975 = mean(Quant_score > 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    qa_table$Sheet <- sheet
    qa_table$ROI <- roi

    output_prefix <- paste0(sheet, "_", roi)
    write.csv(score_output, file.path(sheet_output, paste0(output_prefix, "_new_site_scores.csv")), row.names = FALSE)
    write.csv(calibration_table, file.path(sheet_output, paste0(output_prefix, "_site_calibration.csv")), row.names = FALSE)
    write.csv(qa_table, file.path(sheet_output, paste0(output_prefix, "_site_calibration_QA.csv")), row.names = FALSE)

    calibrated_results <- results
    calibrated_results$site_calibration <- calibration_table
    calibrated_results$site_calibration_QA <- qa_table
    calibrated_results$new_site_score_output <- score_output
    calibrated_results$new_site_scoring_policy <- list(
      calibrated_status = "calibrated_from_healthy_controls",
      insufficient_control_status = "uncalibrated_neutral_fallback",
      minimum_controls = minimum_controls,
      diagnosis_column = diagnosis_column,
      healthy_labels = healthy_labels,
      neutral_site = neutral$neutral_site,
      neutral_site_method = neutral$method,
      percentile_definition = "zzz_cent GG cumulative distribution function",
      z_score_definition = "zzz_cent normalized Z score"
    )
    saveRDS(
      calibrated_results,
      file.path(sheet_output, paste0(output_prefix, "_new_site_calibrated_results.rds"))
    )

    all_calibration_tables[[paste(sheet, roi, sep = "__")]] <- calibration_table
    all_qa_tables[[paste(sheet, roi, sep = "__")]] <- qa_table
  }
}

write.csv(
  bind_rows(all_calibration_tables),
  file.path(output_root, "ALL_ROIS_new_site_calibration_summary.csv"),
  row.names = FALSE
)
write.csv(
  bind_rows(all_qa_tables),
  file.path(output_root, "ALL_ROIS_new_site_calibration_QA.csv"),
  row.names = FALSE
)

message("New-site calibration and scoring completed.")

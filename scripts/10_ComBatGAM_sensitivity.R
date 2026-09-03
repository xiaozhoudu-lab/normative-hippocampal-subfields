# ============================================================
# 10_ComBatGAM_sensitivity.R
# Method S2. Sensitivity analyses for site-related
# harmonization and scale effects
#
# Part 1. HC-only ComBat-GAM harmonization of 24 bilateral
#         hippocampal measures, followed by GAMLSS refitting
#         without site/scanner random effects.
# Part 2. Empirical site-specific scale sensitivity analysis
#         using sigma_adjusted = sigma * site_scale_factor.
# ============================================================

rm(list = ls())

# Usage:
# Rscript 10_ComBatGAM_sensitivity.R \
#   <source_dir> <mr_subfield_workbook.xlsx> <clinical_data.csv> \
#   <primary_model_root> <output_dir> [analysis_part] \
#   [minimum_site_hc] [age_spline_k] [seed]
#
# analysis_part: all (default), combatgam, or site_scale

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 10_ComBatGAM_sensitivity.R <source_dir> ",
    "<mr_subfield_workbook.xlsx> <clinical_data.csv> ",
    "<primary_model_root> <output_dir> [analysis_part] ",
    "[minimum_site_hc] [age_spline_k] [seed]"
  )
}

source_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
mr_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
clinical_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
primary_model_root <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
output_dir <- args[[5]]
analysis_part <- if (length(args) >= 6) tolower(args[[6]]) else "all"
minimum_site_hc <- if (length(args) >= 7) as.integer(args[[7]]) else 10L
age_spline_k <- if (length(args) >= 8) as.integer(args[[8]]) else 10L
random_seed <- if (length(args) >= 9) as.integer(args[[9]]) else 123L

if (!(analysis_part %in% c("all", "combatgam", "site_scale"))) {
  stop("analysis_part must be all, combatgam, or site_scale")
}
if (!is.finite(minimum_site_hc) || minimum_site_hc < 2) {
  stop("minimum_site_hc must be an integer of at least 2")
}
if (!is.finite(age_spline_k) || age_spline_k < 3) {
  stop("age_spline_k must be an integer of at least 3")
}
if (!is.finite(random_seed)) stop("seed must be an integer")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(gamlss)
  library(gamlss.add)
  library(gamlss.dist)
  library(mgcv)
  library(ComBatFamily)
  library(ggplot2)
})

for (file in c(
  "100.common-variables.r", "101.common-functions.r",
  "102.gamlss-recode.r", "300.variables.r",
  "301.functions.r", "ZZZ_function.R"
)) {
  source(file.path(source_dir, file))
}

set.seed(random_seed)

age_min <- 4
age_max <- 85
scale_lower_bound <- 0.5
scale_upper_bound <- 2.0

sheets <- c(
  "lh.hipposubfields.vol.table",
  "rh.hipposubfields.vol.table"
)

rois <- c(
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

normalize_sex <- function(x) {
  raw <- tolower(safe_chr(x))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("female", "f", "woman", "0", "2")] <- "Female"
  out[raw %in% c("male", "m", "man", "1")] <- "Male"
  out
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
  stop("Input requires match_id or matching FreeSurfer Path2/Path3 columns")
}

make_site_id <- function(df) {
  existing <- if ("Site_ZZZ" %in% names(df)) safe_chr(df$Site_ZZZ) else rep("", nrow(df))
  components <- if (all(c("Province", "Center", "Manufacturer") %in% names(df))) {
    paste0(safe_chr(df$Province), safe_chr(df$Center), safe_chr(df$Manufacturer))
  } else {
    rep("", nrow(df))
  }
  ifelse(existing != "", existing, components)
}

feature_key <- function(sheet, roi) {
  prefix <- if (startsWith(sheet, "lh.")) "L" else "R"
  paste0(prefix, "_", roi)
}

side_name <- function(sheet) ifelse(startsWith(sheet, "lh."), "left", "right")

safe_cor <- function(x, y, method = "spearman") {
  x <- safe_num(x)
  y <- safe_num(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

safe_mae <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  mean(abs(x[ok] - y[ok]))
}

safe_rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

winsorize_scale <- function(x, lower = 0.5, upper = 2.0) {
  x <- safe_num(x)
  x[!is.finite(x)] <- 1
  pmin(pmax(x, lower), upper)
}

reset_gamlss_formula_environments <- function(model, env = .GlobalEnv) {
  formula_slots <- names(model)[vapply(model, inherits, logical(1), what = "formula")]
  for (slot in formula_slots) environment(model[[slot]]) <- env
  term_slots <- names(model)[vapply(model, inherits, logical(1), what = "terms")]
  for (slot in term_slots) attr(model[[slot]], ".Environment") <- env
  model
}

prime_gamlss_prediction_environment <- function(model, training_data, m0 = NULL) {
  assign("data1", training_data, envir = .GlobalEnv)
  assign("object", model, envir = .GlobalEnv)
  if (!is.null(m0)) assign("m0", m0, envir = .GlobalEnv)
  for (pkg_env in list(asNamespace("gamlss"), .GlobalEnv)) {
    for (env_name in c("gamlss.environment", "gamlss.env")) {
      if (!exists(env_name, envir = pkg_env, inherits = FALSE)) next
      pred_env <- get(env_name, envir = pkg_env, inherits = FALSE)
      if (!is.environment(pred_env)) next
      try(assign("object", model, envir = pred_env), silent = TRUE)
      try(assign("data1", training_data, envir = pred_env), silent = TRUE)
      if (!is.null(m0)) try(assign("m0", m0, envir = pred_env), silent = TRUE)
    }
  }
  invisible(TRUE)
}

predict_parameters <- function(model, m0, training_data, scoring_data) {
  model <- reset_gamlss_formula_environments(model)
  training_data <- as.data.frame(training_data)
  scoring_data <- as.data.frame(scoring_data)

  training_data$Sex <- droplevels(factor(as.character(training_data$Sex)))
  scoring_data$Sex <- factor(as.character(scoring_data$Sex), levels = levels(training_data$Sex))
  if (anyNA(scoring_data$Sex)) stop("Sex level absent from model training data")

  if ("Site_ZZZ" %in% names(training_data) && "Site_ZZZ" %in% names(scoring_data)) {
    training_data$Site_ZZZ <- droplevels(factor(as.character(training_data$Site_ZZZ)))
    scoring_data$Site_ZZZ <- factor(
      as.character(scoring_data$Site_ZZZ),
      levels = levels(training_data$Site_ZZZ)
    )
  }

  prime_gamlss_prediction_environment(model, training_data, m0)
  list(
    mu = as.numeric(predict(model, newdata = scoring_data, type = "response", what = "mu")),
    sigma = as.numeric(predict(model, newdata = scoring_data, type = "response", what = "sigma")),
    nu = as.numeric(predict(model, newdata = scoring_data, type = "response", what = "nu"))
  )
}

scores_from_parameters <- function(y, mu, sigma, nu) {
  probability <- suppressWarnings(gamlss.dist::pGG(
    q = y, mu = mu, sigma = sigma, nu = nu
  ))
  probability <- pmin(pmax(as.numeric(probability), 1e-12), 1 - 1e-12)
  list(
    Z_score = stats::qnorm(probability),
    Quant_score = probability
  )
}

find_primary_model <- function(sheet, roi) {
  target_name <- paste0(sheet, "_", roi, "_loop_our_model.rds")
  direct <- file.path(primary_model_root, sheet, target_name)
  if (file.exists(direct)) return(direct)
  candidates <- list.files(
    primary_model_root,
    pattern = paste0("^", gsub("\\.", "\\\\.", target_name), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  candidates <- candidates[
    !grepl("ICV|ComBat|Sensitivity|EmpiricalSiteScale|ScaleEffect", candidates, ignore.case = TRUE)
  ]
  if (length(candidates) != 1) {
    stop(
      "Expected one primary model for ", sheet, " / ", roi,
      "; found ", length(candidates)
    )
  }
  candidates[[1]]
}

read_hc_24roi <- function() {
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_clinical <- c("Age", "Sex", "Diagnosis")
  missing_clinical <- setdiff(required_clinical, names(clinical))
  if (length(missing_clinical) > 0) {
    stop("Clinical data missing: ", paste(missing_clinical, collapse = ", "))
  }
  clinical$match_id <- make_match_id(clinical)
  clinical$Site_ZZZ <- make_site_id(clinical)
  clinical$Age <- safe_num(clinical$Age)
  clinical$Sex <- normalize_sex(clinical$Sex)
  clinical$Diagnosis <- toupper(safe_chr(clinical$Diagnosis))
  clinical <- clinical[clinical$match_id != "", , drop = FALSE]
  if (anyDuplicated(clinical$match_id)) stop("Clinical data contain duplicated match_id")

  if ("Database_included" %in% names(clinical)) {
    clinical <- clinical[safe_num(clinical$Database_included) == 1, , drop = FALSE]
  }
  if ("baseline" %in% names(clinical)) {
    clinical <- clinical[!is.na(clinical$baseline), , drop = FALSE]
  }
  if ("Data_baseline" %in% names(clinical)) {
    clinical <- clinical[!is.na(clinical$Data_baseline), , drop = FALSE]
  }

  read_side <- function(sheet, prefix) {
    side_data <- as.data.frame(read_excel(mr_file, sheet = sheet))
    side_data$match_id <- make_match_id(side_data)
    side_data <- side_data[side_data$match_id != "", , drop = FALSE]
    side_data <- side_data[!duplicated(side_data$match_id), , drop = FALSE]
    missing_roi <- setdiff(rois, names(side_data))
    if (length(missing_roi) > 0) {
      stop("MRI sheet ", sheet, " missing: ", paste(missing_roi, collapse = ", "))
    }
    out <- side_data[, c("match_id", rois), drop = FALSE]
    names(out)[-1] <- paste0(prefix, "_", rois)
    out
  }

  left <- read_side(sheets[[1]], "L")
  right <- read_side(sheets[[2]], "R")
  wide <- clinical %>%
    inner_join(left, by = "match_id") %>%
    inner_join(right, by = "match_id")

  feature_columns <- c(paste0("L_", rois), paste0("R_", rois))
  for (feature in feature_columns) wide[[feature]] <- safe_num(wide[[feature]])

  complete_positive <- apply(
    wide[, feature_columns, drop = FALSE],
    1,
    function(x) all(is.finite(x) & x > 0)
  )
  wide <- wide[
    wide$Diagnosis == "HC" &
      is.finite(wide$Age) & wide$Age >= age_min & wide$Age <= age_max &
      wide$Sex %in% c("Female", "Male") &
      wide$Site_ZZZ != "" & complete_positive,
    , drop = FALSE
  ]
  if (nrow(wide) == 0) stop("No complete HC records remain for ComBat-GAM")

  site_counts_before <- as.data.frame(table(wide$Site_ZZZ), stringsAsFactors = FALSE)
  names(site_counts_before) <- c("Site_ZZZ", "n_HC_complete_24ROI")
  eligible_sites <- site_counts_before$Site_ZZZ[
    site_counts_before$n_HC_complete_24ROI >= minimum_site_hc
  ]
  wide <- wide[wide$Site_ZZZ %in% eligible_sites, , drop = FALSE]
  if (nrow(wide) == 0 || length(eligible_sites) < 2) {
    stop("ComBat-GAM requires at least two sites meeting minimum_site_hc")
  }
  wide <- wide[order(wide$Site_ZZZ, wide$Age, wide$Sex, wide$match_id), , drop = FALSE]
  rownames(wide) <- wide$match_id

  site_counts_after <- as.data.frame(table(wide$Site_ZZZ), stringsAsFactors = FALSE)
  names(site_counts_after) <- c("Site_ZZZ", "n_HC_used_by_ComBatGAM")
  list(
    wide = wide,
    feature_columns = feature_columns,
    site_counts_before = site_counts_before,
    site_counts_after = site_counts_after
  )
}

fit_fp_no_site <- function(mu_poly, sigma_poly, input_data, model_control) {
  gamlss(
    formula = as.formula(paste0(
      "feature ~ fp(Age, npoly = ", mu_poly, ") + Sex"
    )),
    sigma.formula = as.formula(paste0(
      "feature ~ fp(Age, npoly = ", sigma_poly, ") + Sex"
    )),
    family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
    control = model_control,
    data = input_data
  )
}

select_fp_no_site <- function(training_data, seed) {
  set.seed(seed)
  child <- training_data[training_data$Age <= 18, , drop = FALSE]
  adult <- training_data[training_data$Age > 18 & training_data$Age < 70, , drop = FALSE]
  older <- training_data[training_data$Age >= 70, , drop = FALSE]
  if (nrow(adult) > 0) {
    adult <- adult[sample(seq_len(nrow(adult)), ceiling(0.30 * nrow(adult))), , drop = FALSE]
  }
  selection_data <- rbind(child, adult, older)
  control <- gamlss.control()

  rows <- list()
  index <- 0L
  for (mu_poly in 1:3) {
    for (sigma_poly in 1:3) {
      index <- index + 1L
      fitted <- try(
        fit_fp_no_site(mu_poly, sigma_poly, selection_data, control),
        silent = TRUE
      )
      rows[[index]] <- data.frame(
        mu_poly = mu_poly,
        sigma_poly = sigma_poly,
        mu_random = 0L,
        sigma_random = 0L,
        Global_Deviance = if (inherits(fitted, "try-error")) NA_real_ else fitted$G.deviance,
        AIC = if (inherits(fitted, "try-error")) NA_real_ else fitted$aic,
        SBC = if (inherits(fitted, "try-error")) NA_real_ else fitted$sbc,
        BIC = if (inherits(fitted, "try-error")) NA_real_ else BIC(fitted),
        status = if (inherits(fitted, "try-error")) "failed" else "OK"
      )
    }
  }
  table <- bind_rows(rows)
  valid <- which(is.finite(table$BIC))
  if (length(valid) == 0) {
    selected_mu <- 2L
    selected_sigma <- 2L
  } else {
    best <- valid[which.min(table$BIC[valid])]
    selected_mu <- table$mu_poly[[best]]
    selected_sigma <- table$sigma_poly[[best]]
  }
  list(
    table = table,
    selected_mu = selected_mu,
    selected_sigma = selected_sigma,
    control = control
  )
}

fit_combatgam_gamlss <- function(training_data, sheet, roi, model_dir) {
  training_data <- training_data[
    is.finite(training_data$feature) & training_data$feature > 0,
    , drop = FALSE
  ]
  mean_feature <- mean(training_data$feature)
  sd_feature <- stats::sd(training_data$feature)
  training_data <- training_data[
    training_data$feature > mean_feature - 3 * sd_feature &
      training_data$feature < mean_feature + 3 * sd_feature,
    , drop = FALSE
  ]
  training_data$Sex <- factor(training_data$Sex, levels = c("Female", "Male"))
  training_data$Site_ZZZ <- factor(training_data$Site_ZZZ)
  training_data$tem_feature <- training_data$feature
  training_data <- training_data[order(training_data$Age), , drop = FALSE]
  if (nrow(training_data) < 50) stop("Too few HC rows after ROI QC")

  selected <- select_fp_no_site(training_data, random_seed)
  m0 <- fit_fp_no_site(
    selected$selected_mu,
    selected$selected_sigma,
    training_data,
    selected$control
  )
  m2 <- gamlss(
    formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex,
    sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex,
    family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
    control = selected$control,
    data = training_data
  )

  parameters <- predict_parameters(m2, m0, training_data, training_data)
  scores <- scores_from_parameters(
    training_data$feature, parameters$mu, parameters$sigma, parameters$nu
  )
  score_output <- data.frame(
    match_id = rownames(training_data),
    Age = training_data$Age,
    Sex = as.character(training_data$Sex),
    Site_ZZZ = as.character(training_data$Site_ZZZ),
    harmonized_volume = training_data$feature,
    mu = parameters$mu,
    sigma = parameters$sigma,
    nu = parameters$nu,
    Z_score = scores$Z_score,
    Quant_score = scores$Quant_score,
    stringsAsFactors = FALSE
  )

  age_grid <- seq(age_min, age_max, length.out = 1000)
  curve_grid <- expand.grid(
    Age = age_grid,
    Sex = c("Female", "Male"),
    stringsAsFactors = FALSE
  )
  curve_grid$Sex <- factor(curve_grid$Sex, levels = levels(training_data$Sex))
  curve_parameters <- predict_parameters(m2, m0, training_data, curve_grid)
  for (probability in c(0.005, 0.025, 0.5, 0.975, 0.995)) {
    label <- paste0("centile_", gsub("\\.", "_", format(probability * 100, trim = TRUE)))
    curve_grid[[label]] <- gamlss.dist::qGG(
      probability,
      mu = curve_parameters$mu,
      sigma = curve_parameters$sigma,
      nu = curve_parameters$nu
    )
  }
  curve_grid$mu <- curve_parameters$mu
  curve_grid$sigma <- curve_parameters$sigma
  curve_grid$nu <- curve_parameters$nu

  model_sheet_dir <- file.path(model_dir, sheet)
  dir.create(model_sheet_dir, recursive = TRUE, showWarnings = FALSE)
  result <- list(
    m0 = m0,
    m2 = m2,
    data1 = training_data,
    score_output = score_output,
    curve_data = curve_grid,
    list_fit = selected$table,
    str = sheet,
    i = roi,
    selected_mu_poly = selected$selected_mu,
    selected_sigma_poly = selected$selected_sigma,
    selected_mu_random = 0L,
    selected_sigma_random = 0L,
    age_range = "4-85",
    ComBatGAM_adjusted = TRUE,
    ComBatGAM_formula = paste0("y ~ s(Age, k = ", age_spline_k, ") + Sex"),
    ComBatGAM_minimum_site_HC = minimum_site_hc
  )
  model_file <- file.path(
    model_sheet_dir,
    paste0(sheet, "_", roi, "_loop_our_model_ComBatGAM.rds")
  )
  saveRDS(result, model_file)
  write.csv(
    selected$table,
    file.path(model_sheet_dir, paste0(sheet, "_", roi, "_model_selection_ComBatGAM.csv")),
    row.names = FALSE
  )
  write.csv(
    curve_grid,
    file.path(model_sheet_dir, paste0(sheet, "_", roi, "_centile_curves_ComBatGAM.csv")),
    row.names = FALSE
  )
  list(result = result, model_file = model_file)
}

run_empirical_site_scale <- function() {
  part_dir <- file.path(output_dir, "Empirical_site_scale")
  subject_dir <- file.path(part_dir, "Subject_Level")
  dir.create(subject_dir, recursive = TRUE, showWarnings = FALSE)

  summary_rows <- list()
  scale_rows <- list()
  failure_rows <- list()
  result_index <- 0L
  scale_index <- 0L
  failure_index <- 0L

  for (sheet in sheets) {
    for (roi in rois) {
      model_file <- try(find_primary_model(sheet, roi), silent = TRUE)
      result <- if (inherits(model_file, "try-error")) model_file else try(readRDS(model_file), silent = TRUE)
      if (inherits(result, "try-error") || is.null(result$m2) ||
          is.null(result$m0) || is.null(result$data1)) {
        failure_index <- failure_index + 1L
        failure_rows[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = "primary model unavailable or incomplete"
        )
        next
      }

      training <- as.data.frame(result$data1)
      required <- c("Age", "Sex", "Site_ZZZ", "feature")
      if (length(setdiff(required, names(training))) > 0) {
        failure_index <- failure_index + 1L
        failure_rows[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = "primary data1 lacks required columns"
        )
        next
      }
      training$Age <- safe_num(training$Age)
      training$feature <- safe_num(training$feature)
      training$Sex <- factor(as.character(training$Sex))
      training$Site_ZZZ <- factor(as.character(training$Site_ZZZ))
      valid <- is.finite(training$Age) & is.finite(training$feature) &
        !is.na(training$Sex) & !is.na(training$Site_ZZZ)
      training <- training[valid, , drop = FALSE]

      parameters <- try(
        predict_parameters(result$m2, result$m0, training, training),
        silent = TRUE
      )
      if (inherits(parameters, "try-error")) {
        failure_index <- failure_index + 1L
        failure_rows[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = "primary parameter prediction failed"
        )
        next
      }
      primary_scores <- scores_from_parameters(
        training$feature, parameters$mu, parameters$sigma, parameters$nu
      )
      subject <- data.frame(
        match_id = rownames(training),
        Age = training$Age,
        Sex = as.character(training$Sex),
        Site_ZZZ = as.character(training$Site_ZZZ),
        observed_volume = training$feature,
        mu = parameters$mu,
        sigma = parameters$sigma,
        nu = parameters$nu,
        primary_z = primary_scores$Z_score,
        primary_centile = primary_scores$Quant_score,
        stringsAsFactors = FALSE
      )

      site_scale <- subject %>%
        group_by(Site_ZZZ) %>%
        summarise(
          site_n = n(),
          site_z_mean = mean(primary_z, na.rm = TRUE),
          site_z_sd_raw = stats::sd(primary_z, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          site_scale_factor_raw = ifelse(
            site_n >= minimum_site_hc & is.finite(site_z_sd_raw) & site_z_sd_raw > 0,
            site_z_sd_raw,
            1
          ),
          site_scale_factor = winsorize_scale(
            site_scale_factor_raw, scale_lower_bound, scale_upper_bound
          ),
          scale_factor_used = site_n >= minimum_site_hc &
            is.finite(site_z_sd_raw) & site_z_sd_raw > 0,
          sheet = sheet,
          side = side_name(sheet),
          ROI = roi
        )

      subject <- left_join(
        subject,
        site_scale[, c("Site_ZZZ", "site_n", "site_scale_factor")],
        by = "Site_ZZZ"
      )
      subject$site_scale_factor[!is.finite(subject$site_scale_factor)] <- 1
      subject$sigma_adjusted <- subject$sigma * subject$site_scale_factor
      adjusted_scores <- scores_from_parameters(
        subject$observed_volume,
        subject$mu,
        subject$sigma_adjusted,
        subject$nu
      )
      subject$adjusted_z <- adjusted_scores$Z_score
      subject$adjusted_centile <- adjusted_scores$Quant_score
      subject$side <- side_name(sheet)
      subject$ROI <- roi

      subject_file <- file.path(
        subject_dir,
        paste0(side_name(sheet), "_", roi, "_empirical_site_scale_scores.csv")
      )
      write.csv(subject, subject_file, row.names = FALSE)

      scale_index <- scale_index + 1L
      scale_rows[[scale_index]] <- site_scale
      used_scales <- site_scale$site_scale_factor[site_scale$scale_factor_used]
      result_index <- result_index + 1L
      summary_rows[[result_index]] <- data.frame(
        sheet = sheet,
        side = side_name(sheet),
        ROI = roi,
        n_HC = nrow(subject),
        n_sites = nrow(site_scale),
        n_sites_scale_estimated = sum(site_scale$scale_factor_used),
        centile_spearman = safe_cor(
          subject$primary_centile, subject$adjusted_centile, "spearman"
        ),
        centile_MAE = safe_mae(subject$primary_centile, subject$adjusted_centile),
        centile_MAE_percentage_points = 100 * safe_mae(
          subject$primary_centile, subject$adjusted_centile
        ),
        primary_prop_below_2.5 = mean(subject$primary_centile < 0.025, na.rm = TRUE),
        adjusted_prop_below_2.5 = mean(subject$adjusted_centile < 0.025, na.rm = TRUE),
        primary_prop_above_97.5 = mean(subject$primary_centile > 0.975, na.rm = TRUE),
        adjusted_prop_above_97.5 = mean(subject$adjusted_centile > 0.975, na.rm = TRUE),
        median_site_scale_factor = if (length(used_scales)) median(used_scales) else NA_real_,
        min_site_scale_factor = if (length(used_scales)) min(used_scales) else NA_real_,
        max_site_scale_factor = if (length(used_scales)) max(used_scales) else NA_real_,
        subject_level_file = subject_file,
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- bind_rows(summary_rows)
  scales <- bind_rows(scale_rows)
  write.csv(
    summary,
    file.path(part_dir, "empirical_site_scale_sensitivity_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    scales,
    file.path(part_dir, "empirical_site_scale_factors_by_site.csv"),
    row.names = FALSE
  )
  if (length(failure_rows)) {
    write.csv(
      bind_rows(failure_rows),
      file.path(part_dir, "empirical_site_scale_failures.csv"),
      row.names = FALSE
    )
  }
  list(summary = summary, scales = scales)
}

run_combatgam <- function() {
  part_dir <- file.path(output_dir, "ComBatGAM")
  model_dir <- file.path(part_dir, "Models")
  subject_dir <- file.path(part_dir, "Subject_Level")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(subject_dir, recursive = TRUE, showWarnings = FALSE)

  input <- read_hc_24roi()
  wide <- input$wide
  feature_columns <- input$feature_columns
  pheno <- data.frame(
    match_id = wide$match_id,
    Age = wide$Age,
    Sex = factor(wide$Sex, levels = c("Female", "Male")),
    Site_ZZZ = factor(wide$Site_ZZZ),
    stringsAsFactors = FALSE
  )
  combat_input <- as.matrix(wide[, feature_columns, drop = FALSE])
  rownames(combat_input) <- wide$match_id
  covariates <- pheno[, c("Age", "Sex"), drop = FALSE]

  set.seed(random_seed)
  message(
    "Running ComBat-GAM: y ~ s(Age, k = ", age_spline_k,
    ") + Sex; batch = Site_ZZZ"
  )
  combat_fit <- ComBatFamily::combat_gam(
    data = combat_input,
    bat = pheno$Site_ZZZ,
    covar = covariates,
    formula = y ~ s(Age, k = age_spline_k) + Sex,
    eb = TRUE
  )
  harmonized_matrix <- combat_fit$dat.combat
  if (!identical(dim(harmonized_matrix), dim(combat_input))) {
    stop("ComBat-GAM output dimensions do not match its input")
  }
  colnames(harmonized_matrix) <- feature_columns
  rownames(harmonized_matrix) <- wide$match_id
  harmonized <- as.data.frame(harmonized_matrix, check.names = FALSE)
  harmonized$match_id <- rownames(harmonized_matrix)

  write.csv(pheno, file.path(part_dir, "combatgam_pheno_HC.csv"), row.names = FALSE)
  write.csv(
    wide[, c("match_id", feature_columns)],
    file.path(part_dir, "HC_original_24ROI_wide.csv"),
    row.names = FALSE
  )
  write.csv(
    harmonized[, c("match_id", feature_columns)],
    file.path(part_dir, "HC_ComBatGAM_24ROI_wide.csv"),
    row.names = FALSE
  )
  write.csv(
    input$site_counts_before,
    file.path(part_dir, "combatgam_site_counts_before_filter.csv"),
    row.names = FALSE
  )
  write.csv(
    input$site_counts_after,
    file.path(part_dir, "combatgam_site_counts_after_filter.csv"),
    row.names = FALSE
  )
  saveRDS(combat_fit, file.path(part_dir, "combatgam_fit_24ROI.rds"))

  summary_rows <- list()
  failures <- list()
  summary_index <- 0L
  failure_index <- 0L

  for (sheet in sheets) {
    for (roi in rois) {
      key <- feature_key(sheet, roi)
      roi_data <- data.frame(
        match_id = wide$match_id,
        Age = wide$Age,
        Sex = wide$Sex,
        Site_ZZZ = wide$Site_ZZZ,
        original_volume = safe_num(wide[[key]]),
        feature = safe_num(harmonized[[key]]),
        stringsAsFactors = FALSE,
        row.names = wide$match_id
      )
      valid_harmonized <- is.finite(roi_data$feature) & roi_data$feature > 0
      training <- roi_data[valid_harmonized, , drop = FALSE]

      fitted <- try(
        fit_combatgam_gamlss(training, sheet, roi, model_dir),
        silent = TRUE
      )
      if (inherits(fitted, "try-error")) {
        failure_index <- failure_index + 1L
        failures[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = as.character(fitted),
          stringsAsFactors = FALSE
        )
        next
      }
      combat_result <- fitted$result
      combat_scores <- combat_result$score_output

      primary_file <- try(find_primary_model(sheet, roi), silent = TRUE)
      primary_result <- if (inherits(primary_file, "try-error")) {
        primary_file
      } else {
        try(readRDS(primary_file), silent = TRUE)
      }
      if (inherits(primary_result, "try-error") || is.null(primary_result$m2) ||
          is.null(primary_result$m0) || is.null(primary_result$data1)) {
        failure_index <- failure_index + 1L
        failures[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = "primary comparison model unavailable",
          stringsAsFactors = FALSE
        )
        next
      }

      primary_training <- as.data.frame(primary_result$data1)
      primary_scoring <- roi_data[combat_scores$match_id, , drop = FALSE]
      primary_scoring$feature <- primary_scoring$original_volume
      primary_scoring$tem_feature <- primary_scoring$original_volume
      parameters <- try(
        predict_parameters(
          primary_result$m2,
          primary_result$m0,
          primary_training,
          primary_scoring
        ),
        silent = TRUE
      )
      if (inherits(parameters, "try-error")) {
        failure_index <- failure_index + 1L
        failures[[failure_index]] <- data.frame(
          sheet = sheet, ROI = roi, reason = "primary score prediction failed",
          stringsAsFactors = FALSE
        )
        next
      }
      primary_scores <- scores_from_parameters(
        primary_scoring$original_volume,
        parameters$mu,
        parameters$sigma,
        parameters$nu
      )

      comparison <- data.frame(
        match_id = combat_scores$match_id,
        Age = combat_scores$Age,
        Sex = combat_scores$Sex,
        Site_ZZZ = combat_scores$Site_ZZZ,
        original_volume = primary_scoring$original_volume,
        harmonized_volume = combat_scores$harmonized_volume,
        primary_Z = primary_scores$Z_score,
        combatgam_Z = combat_scores$Z_score,
        primary_centile = primary_scores$Quant_score,
        combatgam_centile = combat_scores$Quant_score,
        side = side_name(sheet),
        ROI = roi,
        stringsAsFactors = FALSE
      )
      subject_file <- file.path(
        subject_dir,
        paste0(side_name(sheet), "_", roi, "_primary_vs_ComBatGAM.csv")
      )
      write.csv(comparison, subject_file, row.names = FALSE)

      summary_index <- summary_index + 1L
      summary_rows[[summary_index]] <- data.frame(
        sheet = sheet,
        side = side_name(sheet),
        ROI = roi,
        n_HC_ComBatGAM = nrow(combat_scores),
        n_compared = nrow(comparison),
        selected_mu_poly = combat_result$selected_mu_poly,
        selected_sigma_poly = combat_result$selected_sigma_poly,
        selected_mu_random = 0L,
        selected_sigma_random = 0L,
        centile_spearman = safe_cor(
          comparison$primary_centile, comparison$combatgam_centile, "spearman"
        ),
        centile_MAE = safe_mae(
          comparison$primary_centile, comparison$combatgam_centile
        ),
        centile_MAE_percentage_points = 100 * safe_mae(
          comparison$primary_centile, comparison$combatgam_centile
        ),
        z_spearman = safe_cor(comparison$primary_Z, comparison$combatgam_Z, "spearman"),
        z_MAE = safe_mae(comparison$primary_Z, comparison$combatgam_Z),
        combatgam_model_file = fitted$model_file,
        subject_level_file = subject_file,
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- bind_rows(summary_rows)
  write.csv(
    summary,
    file.path(part_dir, "Primary_vs_ComBatGAM_summary.csv"),
    row.names = FALSE
  )
  if (length(failures)) {
    write.csv(bind_rows(failures), file.path(part_dir, "Failures.csv"), row.names = FALSE)
  }

  nonpositive <- data.frame(
    feature = feature_columns,
    n_nonpositive_or_NA_after_ComBatGAM = vapply(
      feature_columns,
      function(feature) sum(!is.finite(harmonized[[feature]]) | harmonized[[feature]] <= 0),
      integer(1)
    )
  )
  write.csv(
    nonpositive,
    file.path(part_dir, "combatgam_nonpositive_counts_by_feature.csv"),
    row.names = FALSE
  )
  qc <- data.frame(
    setting = c(
      "harmonization", "formula", "batch", "age_min", "age_max",
      "minimum_site_hc", "age_spline_k", "n_HC", "n_sites", "n_features",
      "GAMLSS_site_random_effect"
    ),
    value = c(
      "ComBat-GAM", paste0("y ~ s(Age, k = ", age_spline_k, ") + Sex"),
      "Site_ZZZ", age_min, age_max, minimum_site_hc, age_spline_k,
      nrow(wide), length(unique(wide$Site_ZZZ)), length(feature_columns), "none"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(qc, file.path(part_dir, "combatgam_QC_summary.csv"), row.names = FALSE)
  list(summary = summary, qc = qc)
}

combat_result <- NULL
scale_result <- NULL
if (analysis_part %in% c("all", "combatgam")) combat_result <- run_combatgam()
if (analysis_part %in% c("all", "site_scale")) scale_result <- run_empirical_site_scale()

if (!is.null(combat_result) && !is.null(scale_result)) {
  combined <- full_join(
    combat_result$summary,
    scale_result$summary,
    by = c("sheet", "side", "ROI"),
    suffix = c("_ComBatGAM", "_site_scale")
  )
  write.csv(
    combined,
    file.path(output_dir, "Method_S2_agreement_summary.csv"),
    row.names = FALSE
  )

  figure_data <- bind_rows(
    combat_result$summary %>% transmute(
      side, ROI, analysis = "ComBat-GAM",
      spearman = centile_spearman,
      MAE_percentage_points = centile_MAE_percentage_points
    ),
    scale_result$summary %>% transmute(
      side, ROI, analysis = "Empirical site scale",
      spearman = centile_spearman,
      MAE_percentage_points = centile_MAE_percentage_points
    )
  )
  figure_data$feature <- paste0(
    ifelse(figure_data$side == "left", "L.", "R."), figure_data$ROI
  )

  agreement_plot <- ggplot(
    figure_data,
    aes(x = feature, y = spearman, color = analysis, group = analysis)
  ) +
    geom_hline(yintercept = 1, color = "grey60", linewidth = 0.4) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.5) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "Agreement with primary normative centiles",
      x = NULL, y = "Spearman correlation", color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(
    file.path(output_dir, "Method_S2_centile_agreement.png"),
    agreement_plot, width = 13, height = 5, dpi = 300
  )
  ggsave(
    file.path(output_dir, "Method_S2_centile_agreement_editable.pdf"),
    agreement_plot, width = 13, height = 5, device = cairo_pdf
  )
}

write.csv(
  data.frame(
    setting = c(
      "analysis_part", "age_range", "minimum_site_hc", "age_spline_k",
      "ComBatGAM_formula", "ComBatGAM_batch", "ComBatGAM_GAMLSS_site_effect",
      "site_scale_definition", "site_scale_small_site_value", "site_scale_bounds",
      "seed"
    ),
    value = c(
      analysis_part, paste0(age_min, "-", age_max), minimum_site_hc,
      age_spline_k, paste0("y ~ s(Age, k = ", age_spline_k, ") + Sex"),
      "Site_ZZZ", "none", "SD of primary HC z-scores within site", "1",
      paste0(scale_lower_bound, "-", scale_upper_bound), random_seed
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "Method_S2_settings.csv"),
  row.names = FALSE
)

message("Method S2 sensitivity analyses completed: ", analysis_part)

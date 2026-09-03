# ============================================================
# 09_ICV_sensitivity.R
# Hippocampal subfield normative modeling with ICV adjustment
# Age range: 4–85 years
# Original pipeline + all sites included
# Original-style figures + original-style RDS structure
# ============================================================

rm(list = ls())

# Usage:
# Rscript 09_ICV_sensitivity.R <source_dir> \
#   <mr_subfield_workbook.xlsx> <clinical_data.csv> \
#   <icv_workbook.xlsx> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    "Usage: Rscript 09_ICV_sensitivity.R <source_dir> ",
    "<mr_subfield_workbook.xlsx> <clinical_data.csv> ",
    "<icv_workbook.xlsx> <output_dir>"
  )
}

library(readxl)
library(dplyr)
library(gamlss)
library(gamlss.add)
library(ggplot2)
library(reshape2)
library(doParallel)
library(foreach)
library(stringr)

datapath <- paste0(
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE),
  "/"
)
MR_datapath <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
clinical_datapath <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
ICV_datapath <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
savepath <- args[[5]]

setwd(datapath)
if (!dir.exists(savepath)) dir.create(savepath, recursive = TRUE)

source(paste0(datapath, "100.common-variables.r"))
source(paste0(datapath, "101.common-functions.r"))
source(paste0(datapath, "102.gamlss-recode.r"))
source(paste0(datapath, "300.variables.r"))
source(paste0(datapath, "301.functions.r"))
source(paste0(datapath, "ZZZ_function.R"))

var <- c(
  "lh.hipposubfields.vol.table",
  "rh.hipposubfields.vol.table"
)

target_features <- c(
  "subiculum",
  "CA1",
  "presubiculum",
  "molecular_layer_HP",
  "GC.ML.DG",
  "CA3",
  "CA4",
  "Hippocampal_tail",
  "parasubiculum",
  "fimbria",
  "HATA",
  "Whole_hippocampus"
)

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# ============================================================
# Read ICV table
# ============================================================

read_icv_table <- function(ICV_datapath) {
  
  sheet_names <- excel_sheets(ICV_datapath)
  ICV_all <- list()
  
  for (sh in sheet_names) {
    
    tmp <- try(read_excel(ICV_datapath, sheet = sh), silent = TRUE)
    if (inherits(tmp, "try-error")) next
    
    tmp <- as.data.frame(tmp)
    
    if (!all(c("Freesurfer_Path2", "Freesurfer_Path3", "EstimatedTotalIntraCranialVol") %in% colnames(tmp))) {
      next
    }
    
    tmp$Freesurfer_Path2 <- safe_chr(tmp$Freesurfer_Path2)
    tmp$Freesurfer_Path3 <- safe_chr(tmp$Freesurfer_Path3)
    
    tmp$match_id <- paste0(tmp$Freesurfer_Path2, tmp$Freesurfer_Path3)
    tmp$ICV <- safe_num(tmp$EstimatedTotalIntraCranialVol)
    
    tmp <- tmp[
      tmp$match_id != "" &
        !is.na(tmp$ICV) &
        is.finite(tmp$ICV) &
        tmp$ICV > 0,
      c("match_id", "ICV")
    ]
    
    ICV_all[[sh]] <- tmp
  }
  
  ICV_df <- bind_rows(ICV_all)
  
  if (nrow(ICV_df) == 0) {
    stop(
      "No valid ICV records found in the supplied workbook. Required columns: ",
      "Freesurfer_Path2, Freesurfer_Path3, and EstimatedTotalIntraCranialVol."
    )
  }

  conflicting_ids <- ICV_df %>%
    group_by(match_id) %>%
    summarise(n_ICV = n_distinct(round(ICV, 6)), .groups = "drop") %>%
    filter(n_ICV > 1)
  if (nrow(conflicting_ids) > 0) {
    stop(
      "Conflicting ICV values were found for ", nrow(conflicting_ids),
      " match_id(s); resolve duplicate records before modeling. Example: ",
      paste(head(conflicting_ids$match_id, 5), collapse = ", ")
    )
  }
  
  ICV_df <- ICV_df[!duplicated(ICV_df$match_id), ]
  
  cat("ICV records loaded:", nrow(ICV_df), "\n")
  
  return(ICV_df)
}

ICV_df <- read_icv_table(ICV_datapath)

# ============================================================
# Helper functions
# ============================================================

# ICV must be included while selecting the fractional-polynomial degrees and
# site random-effect structure, not only after model selection.  The primary
# analysis helpers intentionally omit ICV, so this sensitivity script uses its
# own otherwise-equivalent fitting helpers.
fit_fp_model_icv <- function(mu_poly, sigma_poly, mu_random, sigma_random,
                             input_data, model_control) {
  mu_formula <- paste0(
    "feature ~ fp(Age, npoly = ", mu_poly, ") + Sex + ICV",
    if (mu_random == 1) " + random(Site_ZZZ)" else ""
  )
  sigma_formula <- paste0(
    "feature ~ fp(Age, npoly = ", sigma_poly, ") + Sex",
    if (sigma_random == 1) " + random(Site_ZZZ)" else ""
  )

  gamlss(
    formula = as.formula(mu_formula),
    sigma.formula = as.formula(sigma_formula),
    family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
    control = model_control,
    data = input_data
  )
}

fit_model_icv <- function(num) {
  pars <- list_par[num, ]
  fitted <- fit_fp_model_icv(
    mu_poly = pars$mu_poly,
    sigma_poly = pars$sigma_poly,
    mu_random = pars$mu_random,
    sigma_random = pars$sigma_random,
    input_data = data1,
    model_control = con
  )

  data.frame(
    mu_poly = pars$mu_poly,
    sigma_poly = pars$sigma_poly,
    mu_random = pars$mu_random,
    sigma_random = pars$sigma_random,
    Global_Deviance = fitted$G.deviance,
    AIC = fitted$aic,
    SBC = fitted$sbc,
    BIC = BIC(fitted),
    Selection = num
  )
}

best_fit_icv <- function(sel_mu_poly, sel_sigma_poly, mu_random, sigma_random) {
  fit_fp_model_icv(
    mu_poly = sel_mu_poly,
    sigma_poly = sel_sigma_poly,
    mu_random = mu_random,
    sigma_random = sigma_random,
    input_data = data1,
    model_control = con
  )
}

make_centile_curve_original <- function(model1, data_ref, sex_values = c("Female", "Male"), num_length = 5000) {
  
  age_grid <- seq(
    min(data_ref$Age, na.rm = TRUE),
    max(data_ref$Age, na.rm = TRUE),
    length.out = num_length
  )
  
  get_site_names <- function(model1, what_name) {
    smo <- model1[[paste0(what_name, ".coefSmo")]]
    if (is.null(smo)) return(NULL)
    
    for (k in seq_along(smo)) {
      cc <- smo[[k]]$coef
      if (!is.null(cc)) {
        nm <- names(cc)
        if (!is.null(nm)) {
          hit <- nm[nm %in% levels(data_ref$Site_ZZZ)]
          if (length(hit) > 0) return(hit)
        }
      }
    }
    
    return(NULL)
  }
  
  average_by_age <- function(x) {
    x <- as.numeric(x)
    
    if (length(x) == num_length) {
      return(x)
    }
    
    if (length(x) %% num_length != 0) {
      stop("Prediction length is not divisible by age-grid length.")
    }
    
    mat <- matrix(x, nrow = num_length)
    rowMeans(mat, na.rm = TRUE)
  }
  
  predict_param <- function(what_name) {
    
    site_names <- get_site_names(model1, what_name)
    
    if (!is.null(site_names)) {
      data3 <- list(
        Age = age_grid,
        Sex = sex_values,
        ICV = median(data_ref$ICV, na.rm = TRUE),
        Site_ZZZ = site_names
      )
      data4 <- do.call(what = expand.grid, args = data3)
      data4$Site_ZZZ <- factor(data4$Site_ZZZ, levels = levels(data_ref$Site_ZZZ))
    } else {
      data3 <- list(
        Age = age_grid,
        Sex = sex_values,
        ICV = median(data_ref$ICV, na.rm = TRUE)
      )
      data4 <- do.call(what = expand.grid, args = data3)
    }
    
    data4$Sex <- factor(data4$Sex, levels = levels(data_ref$Sex))
    
    pred <- predict(
      model1,
      newdata = data4,
      type = "response",
      what = what_name
    )
    
    average_by_age(pred)
  }
  
  mu <- predict_param("mu")
  sigma <- predict_param("sigma")
  nu <- predict_param("nu")
  
  p2 <- zzz_cent(
    obj = model1,
    type = c("centiles"),
    mu = mu,
    sigma = sigma,
    nu = nu,
    cent = c(0.5, 2.5, 50, 97.5, 99.5),
    xname = "Age",
    xvalues = age_grid,
    calibration = FALSE,
    lpar = 3
  )
  
  p2 <- as.data.frame(p2)
  colnames(p2)[1:6] <- c(
    "Age",
    "lower99CI",
    "lower95CI",
    "median",
    "upper95CI",
    "upper99CI"
  )
  
  p2$sigma <- sigma
  
  return(p2)
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
    stop("training_data$Site_ZZZ has no factor levels")
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

make_prediction_data <- function(all_data, model1, training_data, extra_covariates = character()) {
  training_data <- as.data.frame(training_data)
  training_data$Site_ZZZ <- droplevels(factor(as.character(training_data$Site_ZZZ)))
  training_data$Sex <- droplevels(factor(as.character(training_data$Sex)))
  neutral <- choose_neutral_site(model1, training_data)
  preferred <- c("Age", "Sex", "Site_ZZZ", extra_covariates, "tem_feature")
  pred_base <- all_data[, preferred, drop = FALSE]
  rownames(pred_base) <- rownames(all_data)
  original_site <- as.character(pred_base$Site_ZZZ)
  site_unknown <- !(original_site %in% levels(training_data$Site_ZZZ))
  site_used <- original_site
  site_used[site_unknown] <- neutral$neutral_site
  pred_base$Site_ZZZ <- factor(site_used, levels = levels(training_data$Site_ZZZ))
  pred_base$Sex <- factor(as.character(pred_base$Sex), levels = levels(training_data$Sex))
  if (anyNA(pred_base$Sex)) stop("Scoring data contain Sex values absent from model training data")
  pred_base$feature <- pred_base$tem_feature
  keep_cols <- names(pred_base)[names(pred_base) %in% names(training_data)]
  pred_data <- pred_base[, keep_cols, drop = FALSE]
  audit_data <- data.frame(
    Site_ZZZ_original = original_site,
    Site_ZZZ_used_for_prediction = as.character(pred_base$Site_ZZZ),
    Site_unknown_for_model = site_unknown,
    Site_scoring_status = ifelse(site_unknown, "uncalibrated_neutral_fallback", "seen_training_site"),
    neutral_site_used = neutral$neutral_site,
    neutral_site_method = neutral$method,
    stringsAsFactors = FALSE,
    row.names = rownames(all_data)
  )
  list(
    pred_data = pred_data,
    audit_data = audit_data,
    neutral_site = neutral$neutral_site,
    neutral_site_method = neutral$method,
    n_unknown_site = sum(site_unknown, na.rm = TRUE)
  )
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

# ============================================================
# Main loop
# ============================================================

for (sheet in var) {
  
  cat("\n==================================================\n")
  cat("Running sheet:", sheet, "\n")
  cat("==================================================\n")
  
  str <- sheet
  
  outdir <- paste0(savepath, "/", str)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  setwd(datapath)
  MRI <- read_excel(MR_datapath, sheet = sheet)
  MRI <- as.data.frame(MRI)
  
  MRI$Freesurfer_Path2 <- safe_chr(MRI$Freesurfer_Path2)
  MRI$Freesurfer_Path3 <- safe_chr(MRI$Freesurfer_Path3)
  MRI$match_id <- paste0(MRI$Freesurfer_Path2, MRI$Freesurfer_Path3)
  
  MRI <- MRI[MRI$match_id != "", ]
  MRI <- MRI[!duplicated(MRI$match_id), ]
  
  missing_rois <- setdiff(target_features, colnames(MRI))
  if (length(missing_rois) > 0) {
    stop(paste0(
      "These ROIs are not found in sheet ", sheet, ": ",
      paste(missing_rois, collapse = ", ")
    ))
  }
  
  Z_data <- list()
  Quant_data <- list()
  
  for (i in target_features) {
    
    cat("\n--------------------------------------------------\n")
    cat("Running ROI:", i, "| sheet:", sheet, "\n")
    cat("--------------------------------------------------\n")
    
    setwd(datapath)
    data1 <- read.csv(clinical_datapath, header = TRUE, stringsAsFactors = FALSE)
    data1 <- as.data.frame(data1)
    
    data1$Freesufer_Path2 <- safe_chr(data1$Freesufer_Path2)
    data1$Freesufer_Path3 <- safe_chr(data1$Freesufer_Path3)
    data1$match_id <- paste0(data1$Freesufer_Path2, data1$Freesufer_Path3)
    
    data1 <- data1[data1$match_id != "", ]
    data1 <- data1[!duplicated(data1$match_id), ]
    
    if ("Site_ZZZ" %in% colnames(data1)) {
      data1$Site_ZZZ_original <- safe_chr(data1$Site_ZZZ)
    } else {
      data1$Site_ZZZ_original <- ""
    }
    
    data1$Site_ZZZ_from_components <- paste0(
      safe_chr(data1$Province),
      safe_chr(data1$Center),
      safe_chr(data1$Manufacturer)
    )
    
    data1$Site_ZZZ <- ifelse(
      data1$Site_ZZZ_original != "",
      data1$Site_ZZZ_original,
      data1$Site_ZZZ_from_components
    )
    
    data1$Database_included <- 1
    
    matched_merged <- data1 %>%
      inner_join(
        MRI[, c("match_id", "Freesurfer_Path2", "Freesurfer_Path3", i)],
        by = "match_id",
        suffix = c("_clinical", "_MRI")
      ) %>%
      left_join(
        ICV_df,
        by = "match_id"
      )
    
    colnames(matched_merged)[colnames(matched_merged) == i] <- "tem_feature"
    
    matched_merged$Age <- safe_num(matched_merged$Age)
    matched_merged$tem_feature <- safe_num(matched_merged$tem_feature)
    matched_merged$ICV <- safe_num(matched_merged$ICV)
    
    matched_merged <- as.data.frame(matched_merged)
    rownames(matched_merged) <- matched_merged$match_id
    
    setwd(outdir)
    
    write.csv(
      matched_merged,
      file = paste0(str, "_", i, "_01_matched_merged_all_rows_ICV.csv"),
      row.names = FALSE
    )
    
    matched_success <- matched_merged[, c(
      "match_id",
      "Freesufer_Path2",
      "Freesufer_Path3",
      "Diagnosis",
      "Site_ZZZ",
      "ICV"
    )]
    
    write.csv(
      matched_success,
      file = paste0(str, "_", i, "_matched_success_list_ICV.csv"),
      row.names = FALSE
    )
    
    all_data <- matched_merged
    
    if ("baseline" %in% colnames(all_data)) {
      all_data <- all_data[!is.na(all_data$baseline), ]
    }
    
    all_data <- all_data[
      all_data$Database_included == 1 &
        !is.na(all_data$Age) &
        all_data$Age >= 4 &
        all_data$Age <= 85 &
        !is.na(all_data$Sex) &
        all_data$Sex != "" &
        !is.na(all_data$Site_ZZZ) &
        all_data$Site_ZZZ != "" &
        !is.na(all_data$tem_feature) &
        is.finite(all_data$tem_feature) &
        all_data$tem_feature > 0 &
        !is.na(all_data$ICV) &
        is.finite(all_data$ICV) &
        all_data$ICV > 0,
    ]
    
    rownames(all_data) <- all_data$match_id
    all_data_original <- all_data
    
    data1 <- all_data
    
    if ("Data_baseline" %in% colnames(data1)) {
      data1 <- data1[
        !is.na(data1$tem_feature) &
          !is.na(data1$Data_baseline) &
          data1$Diagnosis == "HC" &
          data1$Age >= 4 &
          data1$Age <= 85,
      ]
    } else {
      data1 <- data1[
        !is.na(data1$tem_feature) &
          data1$Diagnosis == "HC" &
          data1$Age >= 4 &
          data1$Age <= 85,
      ]
    }
    
    data1$Site_ZZZ <- as.factor(data1$Site_ZZZ)
    data1$Sex <- as.factor(data1$Sex)
    data1$Sex <- factor(data1$Sex, levels = c("Female", "Male"))
    
    data1 <- data1[order(data1$Age), ]
    data1[, "feature"] <- data1$tem_feature
    all_data[, "feature"] <- all_data$tem_feature
    
    data1 <- data1[!is.na(data1$tem_feature), ]
    data1 <- data1[
      data1$feature > (mean(data1$feature, na.rm = TRUE) - 3 * sd(data1$feature, na.rm = TRUE)) &
        data1$feature < (mean(data1$feature, na.rm = TRUE) + 3 * sd(data1$feature, na.rm = TRUE)),
    ]
    
    data1 <- data1[data1$feature > 0, ]
    data1 <- data1[, c("Age", "Sex", "Site_ZZZ", "ICV", "tem_feature", "feature")]
    data1 <- na.omit(data1)
    
    cat("\nFinal HC training n:", nrow(data1), "\n")
    
    write.csv(
      data1,
      file = paste0(str, "_", i, "_FINAL_HC_training_population_ICV.csv"),
      row.names = TRUE
    )
    
    # ----------------------------
    # Original adult 30% sampling logic
    # ----------------------------
    data1_backup <- data1
    
    set.seed(123)
    data1_child <- data1[data1$Age <= 18, ]
    data1_adult <- data1[data1$Age > 18 & data1$Age < 70, ]
    data1_old <- data1[data1$Age >= 70, ]
    
    if (nrow(data1_adult) > 0) {
      data1_adult_sample <- data1_adult %>% sample_frac(0.3)
    } else {
      data1_adult_sample <- data1_adult
    }
    
    data1 <- rbind(data1_child, data1_adult_sample, data1_old)
    
    list_par <- data.frame(matrix(0, 3 * 3 * 2 * 2, 4))
    colnames(list_par) <- c("mu_poly", "sigma_poly", "mu_random", "sigma_random")
    
    num <- 0
    for (i_poly in 1:3) {
      for (j_poly in 1:3) {
        for (i_rnd_tmp in 0:1) {
          for (j_rnd_tmp in 0:1) {
            num <- num + 1
            list_par[num, 1] <- i_poly
            list_par[num, 2] <- j_poly
            list_par[num, 3] <- i_rnd_tmp
            list_par[num, 4] <- j_rnd_tmp
          }
        }
      }
    }
    
    con <- gamlss.control()
    
    cl <- NULL
    results_try <- try({
      cl <- makeCluster(10)
      registerDoParallel(cl)
      
      my_data <- foreach(
        num = 1:nrow(list_par),
        .combine = rbind,
        .packages = c("gamlss"),
        .export = c(
          "fit_model_icv", "fit_fp_model_icv", "list_par", "data1", "con"
        )
      ) %dopar% fit_model_icv(num)
      
      list_fit <- as.data.frame(my_data)
      print(list_fit)
      
      model_ind <- which.min(list_fit$BIC)
      sel_mu_poly <- list_fit$mu_poly[model_ind]
      sel_sigma_poly <- list_fit$sigma_poly[model_ind]
      i_rnd <- list_fit$mu_random[model_ind]
      j_rnd <- list_fit$sigma_random[model_ind]
    })

    if (!is.null(cl)) try(stopCluster(cl), silent = TRUE)
    
    if (inherits(results_try, "try-error")) {
      sel_mu_poly <- 2
      sel_sigma_poly <- 2
      i_rnd <- 1
      j_rnd <- 1
      con <- gamlss.control(c.crit = 0.01, n.cyc = 2, autostep = FALSE)
      
      list_fit <- data.frame(
        mu_poly = sel_mu_poly,
        sigma_poly = sel_sigma_poly,
        mu_random = i_rnd,
        sigma_random = j_rnd,
        BIC = NA
      )
    }
    
    data1 <- data1_backup
    
    m0 <- best_fit_icv(sel_mu_poly, sel_sigma_poly, i_rnd, j_rnd)
    
    # ----------------------------
    # Final model with ICV in mu formula
    # ----------------------------
    if (i_rnd == 1 & j_rnd == 1) {
      m2 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + ICV + random(Site_ZZZ),
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
      m3 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + ICV + random(Site_ZZZ),
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + random(Site_ZZZ),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
    } else if (i_rnd == 1 & j_rnd == 0) {
      m2 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + ICV + random(Site_ZZZ),
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex,
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
      m3 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + ICV + random(Site_ZZZ),
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
    } else if (i_rnd == 0 & j_rnd == 1) {
      m2 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + ICV,
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
      m3 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + ICV,
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + random(Site_ZZZ),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
    } else {
      m2 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + ICV,
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex,
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
      
      m3 <- gamlss(
        formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + ICV,
        sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)),
        control = con,
        family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
        data = data1
      )
    }
    
    # ----------------------------
    # Generate centile curves
    # ----------------------------
    num_length <- 5000
    
    p2 <- make_centile_curve_original(m3, data1, sex_values = c("Female", "Male"), num_length = num_length)
    p2_all <- make_centile_curve_original(m2, data1, sex_values = c("Female", "Male"), num_length = num_length)
    male_p2 <- make_centile_curve_original(m2, data1, sex_values = c("Male"), num_length = num_length)
    female_p2 <- make_centile_curve_original(m2, data1, sex_values = c("Female"), num_length = num_length)
    
    colnames(p2) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
    colnames(p2_all) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
    colnames(male_p2) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
    colnames(female_p2) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
    
    step_age <- (max(data1$Age) - min(data1$Age)) / num_length
    
    Grad_p2 <- (p2$median[2:nrow(p2)] - p2$median[1:(nrow(p2) - 1)]) / step_age
    p2 <- cbind(p2, data.frame(Gradient1 = c(Grad_p2, Grad_p2[nrow(p2) - 1])))
    
    Female_p2 <- female_p2
    Male_p2 <- male_p2
    
    Grad_Female_p2 <- (Female_p2$median[2:nrow(Female_p2)] - Female_p2$median[1:(nrow(Female_p2) - 1)]) / step_age
    Female_p2 <- cbind(Female_p2, data.frame(Gradient1 = c(Grad_Female_p2, Grad_Female_p2[nrow(Female_p2) - 1])))
    
    Grad_Male_p2 <- (Male_p2$median[2:nrow(Male_p2)] - Male_p2$median[1:(nrow(Male_p2) - 1)]) / step_age
    Male_p2 <- cbind(Male_p2, data.frame(Gradient1 = c(Grad_Male_p2, Grad_Male_p2[nrow(Male_p2) - 1])))
    
    if (!(str_detect(sheet, "thickness"))) {
      scale1 <- 10000
      ylab1 <- "×10^4 mm3"
    } else {
      scale1 <- 1
      ylab1 <- "mm"
    }
    
    # ----------------------------
    # Original-style figures
    # ----------------------------
    plot_x_scale <- scale_x_continuous(
      breaks = c(6, 18, 35, 80),
      labels = c("6 yr", "18 yr", "35 yr", "80 yr")
    )
    
    png(paste0(str, "_", i, "_all_without_sex_stratified_Gradient_ICV.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)
    print(
      ggplot() +
        geom_line(data = p2, aes(x = Age, y = Gradient1 / scale1), color = "#262626", linewidth = 1) +
        labs(title = paste0(i, " ", ylab1, " ICV-adjusted"), x = "", y = "") +
        theme_bw() +
        plot_x_scale
    )
    dev.off()
    
    png(paste0(str, "_", i, "_all_without_sex_stratified_ICV.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)
    print(
      ggplot() +
        geom_point(data = data1[data1$Sex == "Female", ], aes(x = Age, y = tem_feature / scale1),
                   colour = "#E84935", shape = 16, size = 3, alpha = 0.1) +
        geom_point(data = data1[data1$Sex == "Male", ], aes(x = Age, y = tem_feature / scale1),
                   colour = "#4FBBD8", shape = 17, size = 3, alpha = 0.1) +
        geom_line(data = p2, aes(x = Age, y = median / scale1), color = "#262626", linewidth = 1) +
        geom_line(data = p2, aes(x = Age, y = lower99CI / scale1), color = "#262626", linewidth = 1, linetype = "dashed") +
        geom_line(data = p2, aes(x = Age, y = lower95CI / scale1), color = "#262626", linewidth = 1, linetype = "dotted") +
        geom_line(data = p2, aes(x = Age, y = upper95CI / scale1), color = "#262626", linewidth = 1, linetype = "dotted") +
        geom_line(data = p2, aes(x = Age, y = upper99CI / scale1), color = "#262626", linewidth = 1, linetype = "dashed") +
        labs(title = paste0(i, " ", ylab1, " ICV-adjusted"), x = "", y = "") +
        theme_bw() +
        plot_x_scale
    )
    dev.off()
    
    png(paste0(str, "_", i, "_all_with_sex_stratified_Gradient_ICV.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)
    print(
      ggplot() +
        geom_line(data = Female_p2, aes(x = Age, y = Gradient1 / scale1), color = "#E84935", linewidth = 1) +
        geom_line(data = Male_p2, aes(x = Age, y = Gradient1 / scale1), color = "#4FBBD8", linewidth = 1) +
        labs(title = paste0(i, " ", ylab1, " ICV-adjusted"), x = "", y = "") +
        theme_bw() +
        plot_x_scale
    )
    dev.off()
    
    png(paste0(str, "_", i, "_all_with_sex_stratified_sigma_ICV.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)
    print(
      ggplot() +
        geom_line(data = Female_p2, aes(x = Age, y = sigma), color = "#E84935", linewidth = 1) +
        geom_line(data = Male_p2, aes(x = Age, y = sigma), color = "#4FBBD8", linewidth = 1) +
        labs(title = paste0(i, " ", ylab1, " ICV-adjusted"), x = "", y = "") +
        theme_bw() +
        plot_x_scale
    )
    dev.off()
    
    png(paste0(str, "_", i, "_all_with_sex_stratified_ICV.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)
    print(
      ggplot() +
        geom_point(data = data1[data1$Sex == "Female", ], aes(x = Age, y = tem_feature / scale1),
                   colour = "#E84935", shape = 16, size = 3, alpha = 0.1) +
        geom_point(data = data1[data1$Sex == "Male", ], aes(x = Age, y = tem_feature / scale1),
                   colour = "#4FBBD8", shape = 17, size = 3, alpha = 0.1) +
        geom_line(data = Female_p2, aes(x = Age, y = median / scale1), color = "#E84935", linewidth = 1) +
        geom_line(data = Female_p2, aes(x = Age, y = lower95CI / scale1), color = "#E84935", linewidth = 1, linetype = "dotted") +
        geom_line(data = Female_p2, aes(x = Age, y = upper95CI / scale1), color = "#E84935", linewidth = 1, linetype = "dotted") +
        geom_line(data = Male_p2, aes(x = Age, y = median / scale1), color = "#4FBBD8", linewidth = 1) +
        geom_line(data = Male_p2, aes(x = Age, y = lower95CI / scale1), color = "#4FBBD8", linewidth = 1, linetype = "dotted") +
        geom_line(data = Male_p2, aes(x = Age, y = upper95CI / scale1), color = "#4FBBD8", linewidth = 1, linetype = "dotted") +
        labs(title = paste0(i, " ", ylab1, " ICV-adjusted"), x = "", y = "") +
        theme_bw() +
        plot_x_scale
    )
    dev.off()
    
    # ----------------------------
    # Z-score / percentile for all cases
    # ----------------------------
    all_data <- all_data[
      !is.na(all_data$tem_feature) &
        is.finite(all_data$tem_feature) &
        !is.na(all_data$Age) &
        is.finite(all_data$Age) &
        !is.na(all_data$Site_ZZZ) &
        all_data$Site_ZZZ != "" &
        !is.na(all_data$Sex) &
        all_data$Sex != "" &
        !is.na(all_data$ICV) &
        is.finite(all_data$ICV) &
        all_data$ICV > 0,
    ]
    
    model1 <- m2
    pred_pack <- make_prediction_data(all_data, model1, data1, extra_covariates = "ICV")
    all_data1 <- pred_pack$pred_data
    site_audit <- pred_pack$audit_data

    mu <- predict_gamlss_response(model1, all_data1, data1, "mu", m0)
    sigma <- predict_gamlss_response(model1, all_data1, data1, "sigma", m0)
    nu <- predict_gamlss_response(model1, all_data1, data1, "nu", m0)
    
    Z_score_sum <- zzz_cent(
      obj = model1,
      type = c("z-scores"),
      mu = mu,
      sigma = sigma,
      nu = nu,
      xname = "Age",
      xvalues = all_data1$Age,
      yval = all_data1$tem_feature,
      calibration = FALSE,
      lpar = 3
    )
    
    Quant_score_sum <- zzz_cent(
      obj = model1,
      type = c("z-scores"),
      mu = mu,
      sigma = sigma,
      nu = nu,
      xname = "Age",
      xvalues = all_data1$Age,
      yval = all_data1$tem_feature,
      calibration = FALSE,
      lpar = 3,
      cdf = TRUE
    )
    
    Z_score_sum <- data.frame(Z_score_sum)
    colnames(Z_score_sum) <- c("Z_score")
    rownames(Z_score_sum) <- rownames(all_data1)
    
    Quant_score_sum <- data.frame(Quant_score_sum)
    colnames(Quant_score_sum) <- c("Quant_score")
    rownames(Quant_score_sum) <- rownames(all_data1)
    
    Z_data[[i]] <- Z_score_sum
    Quant_data[[i]] <- Quant_score_sum
    
    score_output <- all_data[rownames(all_data1), , drop = FALSE]
    score_output$Z_score <- Z_score_sum$Z_score
    score_output$Quant_score <- Quant_score_sum$Quant_score
    score_output <- cbind(score_output, site_audit[rownames(all_data1), , drop = FALSE])
    
    write.csv(
      score_output,
      file = paste0(str, "_", i, "_all_diagnoses_scores_ICV.csv"),
      row.names = FALSE
    )
    
    peak_table <- data.frame(
      Sheet = sheet,
      ROI = i,
      peakage = p2$Age[which.max(p2$median)],
      peakage_female = Female_p2$Age[which.max(Female_p2$median)],
      peakage_male = Male_p2$Age[which.max(Male_p2$median)],
      HC_training_n = nrow(data1),
      all_scored_n = nrow(score_output),
      unknown_site_n = pred_pack$n_unknown_site,
      neutral_site_used = pred_pack$neutral_site,
      neutral_site_method = pred_pack$neutral_site_method,
      selected_mu_poly = sel_mu_poly,
      selected_sigma_poly = sel_sigma_poly,
      selected_mu_random = i_rnd,
      selected_sigma_random = j_rnd,
      ICV_adjusted = TRUE,
      ICV_reference_for_curve = median(data1$ICV, na.rm = TRUE)
    )
    
    write.csv(
      peak_table,
      file = paste0(str, "_", i, "_peak_age_ICV.csv"),
      row.names = FALSE
    )
    
    results <- list()
    
    results$Female_p2 <- Female_p2
    results$Male_p2 <- Male_p2
    results$p2 <- p2
    results$peakage <- p2$Age[which.max(p2$median)]
    results$p2_all <- p2_all
    results$m2 <- m2
    results$m0 <- m0
    results$m3 <- m3
    results$list_fit <- list_fit
    results$Zscore <- Z_data
    results$Quant_data <- Quant_data
    results$data1 <- data1
    results$all_data <- all_data1
    results$str <- str
    results$i <- i
    results$all_data_original <- all_data_original
    
    results$matched_success <- matched_success
    results$score_output <- score_output
    results$site_scoring <- list(
      unknown_site_n = pred_pack$n_unknown_site,
      neutral_site_used = pred_pack$neutral_site,
      neutral_site_method = pred_pack$neutral_site_method,
      unknown_site_policy = "uncalibrated_neutral_fallback"
    )
    results$peak_table <- peak_table
    results$selected_mu_poly <- sel_mu_poly
    results$selected_sigma_poly <- sel_sigma_poly
    results$selected_mu_random <- i_rnd
    results$selected_sigma_random <- j_rnd
    results$age_range <- "4-85"
    results$ICV_adjusted <- TRUE
    results$ICV_source_file <- ICV_datapath
    results$ICV_reference_for_curve <- median(data1$ICV, na.rm = TRUE)
    
    saveRDS(
      results,
      paste0(str, "_", i, "_loop_our_model_ICV.rds")
    )
    
    cat("\nFinished ROI:", i, "| sheet:", sheet, "| ICV-adjusted | Age 4-85\n")
  }
}

cat("\nAll done. ICV-adjusted normative models finished for age 4–85.\n")

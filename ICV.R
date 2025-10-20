# ==========================================================
# 01_Normative-model-fit-with-ICV.R
# ==========================================================
# Main script: normative model fitting with ICV adjustment
# (Lifespan hippocampal subfield models, age clamped 4–85)
# 
# This script fits normative lifespan models for hippocampal subfields
# using GAMLSS, incorporating ICV as a covariate and random site effects.
# It also generates site-averaged centile trajectories for both sexes.
# ==========================================================

# ---------------------------- #
# 0. Environment setup
# ---------------------------- #
rm(list = ls())

# ======== User-defined paths (modify these before running) ========
datapath       <- "C:/Users/duxiaozhou/Desktop/DATA/CODE/"   # <- Path containing source functions
clinical_datapath <- file.path(datapath, "new_final_list1_update_20250410_update.csv")
MR_datapath      <- file.path(datapath, "MR_Results_subfield.xlsx")
savepath         <- "C:/Users/duxiaozhou/Desktop/DATA/Results/cm1015"  # <- Output directory

if (!dir.exists(savepath)) dir.create(savepath, recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# 1. Parameters
# ---------------------------- #
AGE_MIN <- 4
AGE_MAX <- 85
cat(sprintf("Running normative model fitting for age range %d–%d\n", AGE_MIN, AGE_MAX))

# ---------------------------- #
# 2. Load shared functions
# ---------------------------- #
source(file.path(datapath, "100.common-variables.r"))
source(file.path(datapath, "101.common-functions.r"))
source(file.path(datapath, "ZZZ_function.R"))
source(file.path(datapath, "300.variables.r"))
source(file.path(datapath, "301.functions.r"))

suppressPackageStartupMessages({
  library(readxl); library(stringr); library(reshape2); library(dplyr)
  library(gamlss); library(ggplot2); library(ggpubr)
  library(doParallel); library(foreach)
})

# ---------------------------- #
# 3. Configuration
# ---------------------------- #
# Only left hippocampal subfields for demonstration
var <- c("lh.hipposubfields.vol.table")

# ---------------------------- #
# 4. Load ICV (aseg.vol.table)
# ---------------------------- #
aseg_raw <- read_excel(MR_datapath, sheet = "aseg.vol.table")
aseg_raw <- subset(
  aseg_raw,
  !Freesurfer_Path2 %in% c("Epilepsy_dicom_info_nii", "ET", "Guojibu_HC_nii", "ZUOXINIAN_nii") &
    !is.na(Freesurfer_Path3)
)
rownames(aseg_raw) <- paste0(aseg_raw$Freesurfer_Path2, aseg_raw$Freesurfer_Path3)

icv_col <- intersect(c("EstimatedTotalIntraCranialVol", "eTIV"), colnames(aseg_raw))[1]
if (is.na(icv_col)) stop("❌ ICV column not found (Expected: EstimatedTotalIntraCranialVol or eTIV).")

# ---------------------------- #
# 5. Function: site-averaged centile generation
# ---------------------------- #
make_centiles_avg_sites <- function(model, data_fit, sex_vec = NULL, n = 5000, ICV_fix = NULL,
                                    age_min = AGE_MIN, age_max = AGE_MAX) {
  stopifnot(is.data.frame(data_fit))
  ages <- seq(age_min, age_max, length.out = n)
  data_fit$Sex <- factor(data_fit$Sex, levels = c("Female", "Male"))
  data_fit$Site_ZZZ <- droplevels(factor(trimws(as.character(data_fit$Site_ZZZ))))
  site_levels <- levels(data_fit$Site_ZZZ)
  include_site <- length(site_levels) > 0
  if (is.null(ICV_fix)) ICV_fix <- median(data_fit$ICV, na.rm = TRUE)
  
  # Grid construction
  build_grid <- function(sex_vec_local) {
    if (is.null(sex_vec_local)) {
      expand.grid(Age = ages, Site_ZZZ = site_levels, ICV = ICV_fix, KEEP.OUT.ATTRS = FALSE)
    } else {
      expand.grid(
        Age = ages,
        Sex = factor(sex_vec_local, levels = c("Female", "Male")),
        Site_ZZZ = site_levels,
        ICV = ICV_fix,
        KEEP.OUT.ATTRS = FALSE
      )
    }
  }
  nd_mu <- build_grid(sex_vec)
  mu_vec <- as.numeric(predict(model, newdata = nd_mu, type = "response", what = "mu"))
  sigma_vec <- as.numeric(predict(model, newdata = nd_mu, type = "response", what = "sigma"))
  nu_vec <- tryCatch(as.numeric(predict(model, newdata = nd_mu, type = "response", what = "nu")),
                     error = function(e) rep(NA_real_, nrow(nd_mu)))
  
  seg <- length(mu_vec) / length(ages)
  avg_over_segments <- function(v) {
    if (seg <= 1) return(v)
    out <- v[1:length(ages)]
    for (k in 2:seg) out <- out + v[((k-1)*length(ages)+1):(k*length(ages))]
    out / seg
  }
  mu <- avg_over_segments(mu_vec)
  sigma <- avg_over_segments(sigma_vec)
  nu <- avg_over_segments(nu_vec)
  
  p <- zzz_cent(
    obj = model, type = "centiles",
    mu = mu, sigma = sigma, nu = nu,
    cent = c(0.5, 2.5, 50, 97.5, 99.5),
    xname = "Age", xvalues = ages, calibration = FALSE, lpar = 3
  )
  p[, "sigma"] <- sigma
  colnames(p) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
  return(p)
}

# ---------------------------- #
# 6. Main modeling loop
# ---------------------------- #
for (sheet in var) {
  cat(sprintf("Processing sheet: %s\n", sheet))
  MRI <- read_excel(MR_datapath, sheet = sheet)
  MRI <- subset(MRI,
                !Freesurfer_Path2 %in% c("Epilepsy_dicom_info_nii", "ET", "Guojibu_HC_nii", "ZUOXINIAN_nii") &
                  !is.na(Freesurfer_Path3))
  rownames(MRI) <- paste0(MRI$Freesurfer_Path2, MRI$Freesurfer_Path3)
  
  if (!str_detect(sheet, "hipposubfields"))
    stop("This script is intended only for *hipposubfields* tables.")
  
  tem_feature <- colnames(MRI)[2:13]
  
  # Load clinical data
  data_clin <- read.csv(clinical_datapath)
  data_clin$Site_ZZZ <- paste0(data_clin$Province, data_clin$Center, data_clin$Manufacturer)
  data_clin$Euler <- rowSums(data_clin[, c("euler_number_l", "euler_number_r")], na.rm = TRUE)
  rownames(data_clin) <- paste0(data_clin$Freesufer_Path2, data_clin$Freesufer_Path3)
  
  # Loop over subfields
  for (feat in tem_feature) {
    prefix <- ifelse(sheet == "lh.hipposubfields.vol.table", "lh_", "rh_")
    rds_file <- file.path(savepath, paste0(prefix, feat, "_loop_our_model.rds"))
    if (file.exists(rds_file)) {
      cat(sprintf("→ Skipping existing model: %s\n", rds_file)); next
    }
    
    inter_row <- intersect(rownames(data_clin), rownames(MRI))
    if (length(inter_row) == 0) {
      warning("No matched subjects, skipping: ", feat); next
    }
    
    data1 <- cbind(data_clin[inter_row, ], MRI[inter_row, feat, drop = FALSE])
    colnames(data1)[ncol(data1)] <- "tem_feature"
    data1$feature <- data1$tem_feature
    data1$ICV <- as.numeric(aseg_raw[inter_row, icv_col])
    data1$Site_ZZZ <- factor(trimws(data1$Site_ZZZ))
    data1$Sex <- factor(data1$Sex, levels = c("Female", "Male"))
    data1 <- data1[order(data1$Age), ]
    
    # Remove outliers (±3 SD)
    mu_f <- mean(data1$feature, na.rm = TRUE)
    sd_f <- sd(data1$feature, na.rm = TRUE)
    data1 <- subset(data1, feature > mu_f - 3*sd_f & feature < mu_f + 3*sd_f)
    
    # Filter by age
    data1 <- subset(data1, Age >= AGE_MIN & Age <= AGE_MAX)
    if (nrow(data1) < 50) { warning("Too few samples after filtering, skip: ", feat); next }
    
    # Fit models (same logic as original)
    # Best-fit power & random effects determined via parallel search
    # ...
    # [Retain your existing best_fit / fit_model functions here]
    
    # Fit final GAMLSS models (m2: with Sex; m3: without Sex)
    # Compute centiles and save outputs
    # ...
    
    # Save results
    saveRDS(results, rds_file)
    cat(sprintf("✅ Model saved: %s\n", rds_file))
  }
}

cat("✅ Normative model fitting completed successfully.\n")

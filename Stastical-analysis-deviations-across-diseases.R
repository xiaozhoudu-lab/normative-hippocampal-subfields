# ============================================================
# Script: 02_Group_Comparison_PermutationTest.R
# Author: Xiaozhou Du (Capital Medical University)
# Date: 2025-10-20
# Version: 1.1
# Description:
#   Group-wise comparison between patients and healthy controls
#   using stratified Monte Carlo permutation tests (Welch t-statistic)
#   under site × sex × 5-year age band stratification.
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(effectsize)
})

set.seed(2025)

#### ------------------------------------------------------
#### User-configurable paths
#### ------------------------------------------------------
# Input: Excel file with columns described above.
input_path  <- "data/all.xlsx"

# Output: CSV summary of all pairwise HC-vs-disease comparisons.
output_path <- "results/permutation_SexAgeSite_FDR.csv"

# Make sure output dir exists
dir.create("results", showWarnings = FALSE, recursive = TRUE)

#### ------------------------------------------------------
#### Analysis configuration
#### ------------------------------------------------------

# Diseases to test (NMOSD intentionally excluded for public release)
diseases <- c("MCI","AD","PD","CSVD","MS")

# ROI columns to analyze
roi_vars <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)

B <- 5000            # number of permutations
HC_case_cap <- 5     # cap on HC:case ratio within each stratum

#### ------------------------------------------------------
#### Helper functions
#### ------------------------------------------------------

# Format median [IQR] as a string like "12.34 [10.00–15.25]"
median_iqr_str <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (length(x) < 1) return(NA_character_)
  med <- median(x)
  q   <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"),
          med, q[1], q[2])
}

# Permute group labels WITHIN strata levels (stratified permutation)
permute_labels_in_strata <- function(group, strata) {
  new <- group
  for (s in levels(strata)) {
    idx <- which(strata == s)
    if (length(idx) > 1) {
      new[idx] <- sample(group[idx], length(idx), replace = FALSE)
    }
  }
  factor(new, levels = levels(group))
}

# Welch-type t statistic for two-group comparison
welch_t_stat <- function(y, g) {
  g <- droplevels(g)
  if (nlevels(g) != 2) return(NA_real_)
  x <- y[g == levels(g)[1]]
  z <- y[g == levels(g)[2]]
  if (length(x) < 2 || length(z) < 2) return(NA_real_)
  sx <- var(x); sz <- var(z)
  mx <- mean(x); mz <- mean(z)
  nx <- length(x); nz <- length(z)
  if (!is.finite(sx) || !is.finite(sz)) return(NA_real_)
  (mz - mx) / sqrt(sx/nx + sz/nz)
}

#### ------------------------------------------------------
#### Load data
#### ------------------------------------------------------

df <- read_excel(input_path) %>%
  mutate(
    # Harmonize diagnosis labels
    Diagnosis = dplyr::recode(
      Diagnosis,
      "MCI"="MCI",
      "AD"="AD",
      "PD"="PD",
      "SVD"="CSVD",
      "CSVD"="CSVD",
      "MS"="MS",
      "HC"="HC",
      "AQP4Pos_NMOSD"="AQP4Pos_NMOSD"  # will not be analyzed here
    ),
    Sex  = as.factor(Sex),
    Site = as.factor(Site)
  ) %>%
  droplevels()

# Basic sanity check for required columns
required_cols <- c("Age","Sex","Diagnosis","Site")
missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop("Missing required columns in input data: ",
       paste(missing_required, collapse = ", "))
}

#### ------------------------------------------------------
#### Main loop over diseases × ROIs
#### ------------------------------------------------------

all_out <- list()

# Define 5-year age bands for stratification
age_min    <- floor(min(df$Age, na.rm = TRUE) / 5) * 5
age_max    <- ceiling(max(df$Age, na.rm = TRUE) / 5) * 5
age_breaks <- seq(age_min, age_max, by = 5)

for (dis in diseases) {
  for (roi in roi_vars) {
    if (!roi %in% names(df)) next

    # Subset to HC and this disease only
    sub <- df %>%
      dplyr::select(Age, Sex, Diagnosis, Site, !!sym(roi)) %>%
      dplyr::rename(ROI = !!sym(roi)) %>%
      dplyr::filter(Diagnosis %in% c("HC", dis)) %>%
      tidyr::drop_na(Age, Sex, Diagnosis, Site, ROI) %>%
      dplyr::mutate(
        AgeBand = cut(
          Age,
          breaks = age_breaks,
          include.lowest = TRUE,
          right = FALSE
        )
      )

    # Keep only strata (Sex × AgeBand × Site) where both HC and Disease exist
    sub <- sub %>%
      dplyr::group_by(Sex, AgeBand, Site) %>%
      dplyr::mutate(
        n_HC  = sum(Diagnosis == "HC"),
        n_DIS = sum(Diagnosis == dis)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(n_HC > 0 & n_DIS > 0)

    # If too few samples remain, skip this ROI for this disease
    if (nrow(sub) < 10) next

    # Downsample HC within each (Sex × AgeBand × Site) stratum
    sub <- sub %>%
      dplyr::group_by(Sex, AgeBand, Site) %>%
      dplyr::group_modify(function(d, k) {
        n_dis <- sum(d$Diagnosis == dis)
        cap   <- HC_case_cap * max(1, n_dis)
        hc_id <- which(d$Diagnosis == "HC")
        if (length(hc_id) > cap) {
          keep_hc <- sample(hc_id, cap)
          d <- d[c(keep_hc, which(d$Diagnosis == dis)), , drop = FALSE]
        }
        d
      }) %>%
      dplyr::ungroup()

    # Create strata label and final 2-level group factor
    sub$strata <- interaction(sub$Sex, sub$AgeBand, sub$Site, drop = TRUE)
    sub$group  <- factor(
      ifelse(sub$Diagnosis == "HC", "HC", "Disease"),
      levels = c("HC","Disease")
    )

    # Observed Welch-type t
    t_obs <- welch_t_stat(sub$ROI, sub$group)
    if (!is.finite(t_obs)) next

    # Stratified permutation to get empirical p-value
    t_perm <- numeric(B)
    for (b in seq_len(B)) {
      g_star     <- permute_labels_in_strata(sub$group, sub$strata)
      t_perm[b]  <- welch_t_stat(sub$ROI, g_star)
    }
    p_perm <- (sum(abs(t_perm) >= abs(t_obs)) + 1) / (B + 1)

    # Descriptive stats
    med_iqr_HC  <- median_iqr_str(sub$ROI[sub$group == "HC"])
    med_iqr_DIS <- median_iqr_str(sub$ROI[sub$group == "Disease"])

    # Effect size: Cohen's d (pooled SD) with 95% CI
    d_eff <- tryCatch(
      effectsize::cohens_d(
        ROI ~ group,
        data       = sub,
        pooled_sd  = TRUE,
        ci         = 0.95
      ),
      error = function(e) NULL
    )

    d_val <- if (!is.null(d_eff)) d_eff$Cohens_d[1] else NA_real_
    d_lo  <- if (!is.null(d_eff)) d_eff$CI_low[1]   else NA_real_
    d_hi  <- if (!is.null(d_eff)) d_eff$CI_high[1]  else NA_real_

    # Collect row
    all_out[[length(all_out)+1]] <- data.frame(
      disease        = dis,
      roi            = roi,
      n_HC           = sum(sub$group == "HC"),
      n_Disease      = sum(sub$group == "Disease"),
      median_iqr_HC  = med_iqr_HC,
      median_iqr_dis = med_iqr_DIS,
      t_obs          = t_obs,
      p_perm         = p_perm,
      cohens_d       = d_val,
      d_CI_low       = d_lo,
      d_CI_high      = d_hi,
      stringsAsFactors = FALSE
    )
  }
}

res <- dplyr::bind_rows(all_out)

#### ------------------------------------------------------
#### FDR correction within disease
#### ------------------------------------------------------

if (nrow(res) > 0) {
  res <- res %>%
    dplyr::group_by(disease) %>%
    dplyr::mutate(
      p_fdr = p.adjust(p_perm, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(disease, p_fdr, dplyr::desc(abs(cohens_d)))

  utils::write.csv(res, output_path, row.names = FALSE)
  message("✅ Results saved to: ", output_path)

} else {
  warning(
    "No valid results. Possible reasons:\n",
    "- No overlapping strata (Sex × AgeBand × Site) containing both HC and a given disease;\n",
    "- ROI column names in `roi_vars` not found in input data."
  )
}
# ============================================================
# End of Script
# ============================================================


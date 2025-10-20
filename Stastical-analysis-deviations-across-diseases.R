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

# ---------------------------
# 1. Environment Setup
# ---------------------------
rm(list = ls())
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(effectsize)
})

set.seed(2025)

# ---------------------------
# 2. User Configurations
# ---------------------------
input_path  <- "data/all.xlsx"           # Input Excel file
output_path <- "results/permutation_FDR.csv"  # Output CSV

diseases <- c("MCI","AD","PD","CSVD","MS","AQP4Pos_NMOSD")

roi_vars <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)

B <- 5000        # Number of permutations
HC_case_cap <- 5 # Max HC:Case ratio per stratum

# ---------------------------
# 3. Utility Functions
# ---------------------------

# Median [IQR] string
median_iqr_str <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (length(x) < 1) return(NA_character_)
  med <- median(x)
  q   <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"), med, q[1], q[2])
}

# Permute labels within strata
permute_labels_in_strata <- function(group, strata) {
  new <- group
  for (s in levels(strata)) {
    idx <- which(strata == s)
    if (length(idx) > 1) new[idx] <- sample(group[idx], length(idx))
  }
  factor(new, levels = levels(group))
}

# Welch t-statistic
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

# ---------------------------
# 4. Data Loading
# ---------------------------
message("Reading input data: ", input_path)

df <- read_excel(input_path) %>%
  mutate(
    Diagnosis = recode(Diagnosis,
                       "SVD" = "CSVD",
                       .default = Diagnosis
    ),
    Sex = as.factor(Sex),
    Site_ZZZ = as.factor(Site_ZZZ)
  ) %>%
  droplevels()

stopifnot(all(c("Age","Sex","Site_ZZZ","Diagnosis") %in% names(df)))

# ---------------------------
# 5. Main Analysis Loop
# ---------------------------
all_out <- list()

for (dis in diseases) {
  message("Processing disease: ", dis)
  
  # Define disease-specific age range
  age_min_dis <- floor(min(df$Age[df$Diagnosis == dis], na.rm = TRUE) / 5) * 5
  age_max_dis <- ceiling(max(df$Age[df$Diagnosis == dis], na.rm = TRUE) / 5) * 5
  age_breaks  <- seq(age_min_dis, age_max_dis, by = 5)
  
  for (roi in roi_vars) {
    if (!roi %in% names(df)) next
    
    sub <- df %>%
      select(Age, Sex, Site_ZZZ, Diagnosis, !!sym(roi)) %>%
      rename(ROI = !!sym(roi)) %>%
      filter(Diagnosis %in% c("HC", dis)) %>%
      filter(!(Diagnosis == "HC" & (Age < age_min_dis | Age >= age_max_dis))) %>%
      drop_na()
    
    sub <- sub %>%
      mutate(AgeBand = cut(Age, breaks = age_breaks, include.lowest = TRUE, right = FALSE)) %>%
      group_by(Site_ZZZ, Sex, AgeBand) %>%
      mutate(n_HC = sum(Diagnosis == "HC"),
             n_DIS = sum(Diagnosis == dis)) %>%
      ungroup() %>%
      filter(n_HC > 0 & n_DIS > 0)
    
    if (nrow(sub) < 10) next
    
    # Downsample HC if too many
    sub <- sub %>%
      group_by(Site_ZZZ, Sex, AgeBand) %>%
      group_modify(function(d, k) {
        n_dis <- sum(d$Diagnosis == dis)
        cap   <- HC_case_cap * max(1, n_dis)
        hc_id <- which(d$Diagnosis == "HC")
        if (length(hc_id) > cap) {
          keep_hc <- sample(hc_id, cap)
          d <- d[c(keep_hc, which(d$Diagnosis == dis)), , drop = FALSE]
        }
        d
      }) %>%
      ungroup()
    
    sub$strata <- interaction(sub$Site_ZZZ, sub$Sex, sub$AgeBand, drop = TRUE)
    sub$group  <- factor(ifelse(sub$Diagnosis == "HC", "HC", "Disease"), levels = c("HC","Disease"))
    
    t_obs <- welch_t_stat(sub$ROI, sub$group)
    if (!is.finite(t_obs)) next
    
    t_perm <- replicate(B, {
      g_star <- permute_labels_in_strata(sub$group, sub$strata)
      welch_t_stat(sub$ROI, g_star)
    })
    p_perm <- (sum(abs(t_perm) >= abs(t_obs)) + 1) / (B + 1)
    
    med_iqr_HC  <- median_iqr_str(sub$ROI[sub$group == "HC"])
    med_iqr_DIS <- median_iqr_str(sub$ROI[sub$group == "Disease"])
    
    d_eff <- tryCatch(
      effectsize::hedges_g(ROI ~ group, data = sub, ci = 0.95),
      error = function(e) NULL
    )
    d_val <- if (!is.null(d_eff)) d_eff$Hedges_g[1] else NA_real_
    d_lo  <- if (!is.null(d_eff)) d_eff$CI_low[1]   else NA_real_
    d_hi  <- if (!is.null(d_eff)) d_eff$CI_high[1]  else NA_real_
    
    all_out[[length(all_out)+1]] <- data.frame(
      disease       = dis,
      roi           = roi,
      n_HC          = sum(sub$group == "HC"),
      n_Disease     = sum(sub$group == "Disease"),
      median_iqr_HC = med_iqr_HC,
      median_iqr_dis= med_iqr_DIS,
      t_obs         = t_obs,
      p_perm        = p_perm,
      hedges_g      = d_val,
      g_CI_low      = d_lo,
      g_CI_high     = d_hi,
      stringsAsFactors = FALSE
    )
  }
}

# ---------------------------
# 6. Multiple Comparisons & Output
# ---------------------------
res <- bind_rows(all_out)

if (nrow(res) > 0) {
  res <- res %>%
    group_by(disease) %>%
    mutate(p_fdr = p.adjust(p_perm, method = "BH")) %>%
    ungroup() %>%
    arrange(disease, p_fdr, desc(abs(hedges_g)))
  
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  write.csv(res, output_path, row.names = FALSE)
  message("✅ Results saved to: ", output_path)
} else {
  warning("⚠️ No valid results generated — possibly insufficient overlapping strata.")
}

# ============================================================
# End of Script
# ============================================================

# ==========================================================
# 05_Pointwise-consistency-evaluation.R
# ==========================================================
# Compare Method A vs Method B (no interpolation, point-to-point)
#   Metrics computed:
#     - Pearson correlation (r, p)
#     - MAPE_refB = mean(|A - B| / |B|)
#     - ICC (2,1) with p-value from psych::ICC
#     - CCC (Lin’s concordance correlation coefficient, 95% CI)
#     - Δpeak age difference between curves (years)
#
# Output:
#   - CSV table of all metrics
#   - Boxplot summary (PDF)
#
# Author: Du Xiaozhou
# Version: v1.0 (for GitHub release)
# ==========================================================

rm(list = ls()); options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(readr); library(ggplot2); library(tidyr)
  library(purrr)
  library(psych)   # ICC
  library(epiR)    # CCC
})

# ----------------------------
# 0. USER CONFIGURATION
# ----------------------------
method_A_name <- "NO"      # Label for Method A
method_B_name <- "SITE"    # Label for Method B

# Directories containing .rds files for each method
method_A_dir <- "C:/Users/duxiaozhou/Desktop/DATA/Results/cm_yuan"
method_B_dir <- "C:/Users/duxiaozhou/Desktop/DATA/Results/cm1015"

# Output directory
out_dir <- "C:/Users/duxiaozhou/Desktop/DATA/Results/pointwise_comparison"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Whether to also compute Female / Male strata
compute_sex_strata <- FALSE

# Optional scaling (e.g., divide by 1e4 to match units)
scale_A_by_1e4 <- FALSE
scale_B_by_1e4 <- FALSE


# ----------------------------
# 1. HELPER FUNCTIONS
# ----------------------------

# Extract ROI key from filename
# "lh_CA1_loop_our_model.rds" -> "lh_CA1"
get_roi_key <- function(fname) {
  fn <- basename(fname)
  fn <- sub("\\.rds$", "", fn, ignore.case = TRUE)
  fn <- sub("_loop_our_model$", "", fn, ignore.case = TRUE)
  fn <- sub("_loop_bc_model_bc_peak_age$", "", fn, ignore.case = TRUE)
  fn
}

# Extract curve (overall/Female/Male) from RDS
extract_curve_df <- function(rds_path, split = c("overall", "Female", "Male")) {
  split <- match.arg(split)
  x <- readRDS(rds_path)
  df <- switch(split,
               overall = x$p2,
               Female  = x$Female_p2,
               Male    = x$Male_p2)
  if (is.null(df)) stop("RDS missing required curve: ", rds_path)
  stopifnot(all(c("Age", "median") %in% colnames(df)))
  df %>%
    mutate(Age = as.numeric(Age), median = as.numeric(median)) %>%
    filter(is.finite(Age), is.finite(median)) %>%
    arrange(Age) %>%
    distinct(Age, .keep_all = TRUE)
}

# Align two curves by shared ages
align_pointwise <- function(dfA, dfB) {
  common_age <- intersect(dfA$Age, dfB$Age)
  if (length(common_age) < 3)
    stop("Too few shared ages between A and B.")
  dfA2 <- dfA[dfA$Age %in% common_age, , drop = FALSE][order(dfA$Age[dfA$Age %in% common_age]), ]
  dfB2 <- dfB[dfB$Age %in% common_age, , drop = FALSE][order(dfB$Age[dfB$Age %in% common_age]), ]
  stopifnot(all(dfA2$Age == dfB2$Age))
  list(age = dfA2$Age, A = dfA2$median, B = dfB2$median)
}

# Peak age: prefer RDS$peakage if available
get_peak_age <- function(rds_path, split = c("overall", "Female", "Male")) {
  split <- match.arg(split)
  x <- readRDS(rds_path)
  if (split == "overall" && !is.null(x$peakage)) {
    return(as.numeric(x$peakage))
  } else {
    df <- extract_curve_df(rds_path, split = split)
    if (!nrow(df)) return(NA_real_)
    return(df$Age[which.max(df$median)])
  }
}

# Compute all metrics (Pearson, MAPE, ICC, CCC)
metrics_pointwise <- function(A, B) {
  if (scale_A_by_1e4) A <- A / 10000
  if (scale_B_by_1e4) B <- B / 10000
  
  # Pearson correlation
  ct <- suppressWarnings(cor.test(A, B))
  pearson_r <- unname(ct$estimate)
  pearson_p <- unname(ct$p.value)
  
  # Mean absolute percentage error (ref = B)
  B_safe <- ifelse(B == 0, .Machine$double.eps, B)
  MAPE_refB <- mean(abs(A - B) / abs(B_safe), na.rm = TRUE)
  
  # ICC(2,1)
  ratings <- cbind(A, B)
  icc_out <- try(psych::ICC(ratings, lmer = TRUE), silent = TRUE)
  icc_val <- icc_p <- NA_real_
  if (!inherits(icc_out, "try-error") && "results" %in% names(icc_out)) {
    icc_val <- suppressWarnings(as.numeric(icc_out$results$ICC[2]))
    icc_p   <- suppressWarnings(as.numeric(icc_out$results$p[2]))
  }
  
  # CCC (Lin’s concordance)
  ccc <- try(epiR::epi.ccc(A, B), silent = TRUE)
  ccc_est <- ccc_lcl <- ccc_ucl <- NA_real_
  if (!inherits(ccc, "try-error")) {
    ccc_est <- as.numeric(ccc$rho.c$est)
    ccc_lcl <- as.numeric(ccc$rho.c$lcl)
    ccc_ucl <- as.numeric(ccc$rho.c$ucl)
  }
  
  c(pearson_r = pearson_r, pearson_p = pearson_p,
    MAPE_refB = MAPE_refB,
    ICC = icc_val, ICC_p = icc_p,
    CCC = ccc_est, CCC_lcl = ccc_lcl, CCC_ucl = ccc_ucl)
}

# ----------------------------
# 2. SCAN & MATCH FILES
# ----------------------------
list_rds <- function(dirpath) {
  fs <- list.files(dirpath, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  tibble(file = fs, key = vapply(fs, get_roi_key, FUN.VALUE = character(1)))
}
A_files <- list_rds(method_A_dir)
B_files <- list_rds(method_B_dir)
pairs <- inner_join(A_files, B_files, by = "key", suffix = c("_A", "_B"))
if (nrow(pairs) == 0) stop("No matched ROI RDS files found between directories.")

pairs <- pairs %>%
  mutate(Hemisphere = case_when(
    str_detect(key, "^lh_") ~ "LH",
    str_detect(key, "^rh_") ~ "RH",
    TRUE ~ "NA"
  ),
  ROI = str_replace(key, "^(lh_|rh_)", ""))

# ----------------------------
# 3. MAIN LOOP
# ----------------------------
do_one_split <- function(split = c("overall", "Female", "Male")) {
  split <- match.arg(split)
  pmap_dfr(pairs, function(file_A, key, file_B, Hemisphere, ROI) {
    dfA <- extract_curve_df(file_A, split = split)
    dfB <- extract_curve_df(file_B, split = split)
    aligned <- try(align_pointwise(dfA, dfB), silent = TRUE)
    if (inherits(aligned, "try-error")) {
      return(tibble(Hemisphere, ROI, method_A = method_A_name, method_B = method_B_name,
                    split, pearson_r = NA_real_, pearson_p = NA_real_, MAPE_refB = NA_real_,
                    ICC = NA_real_, ICC_p = NA_real_, CCC = NA_real_,
                    CCC_lcl = NA_real_, CCC_ucl = NA_real_, d_peak_age = NA_real_, n_points = NA_integer_))
    }
    
    mets <- metrics_pointwise(aligned$A, aligned$B)
    pkA <- get_peak_age(file_A, split = split)
    pkB <- get_peak_age(file_B, split = split)
    dpk <- if (is.finite(pkA) & is.finite(pkB)) abs(pkA - pkB) else NA_real_
    
    tibble(Hemisphere, ROI,
           method_A = method_A_name, method_B = method_B_name, split,
           pearson_r = mets["pearson_r"], pearson_p = mets["pearson_p"],
           MAPE_refB = mets["MAPE_refB"],
           ICC = mets["ICC"], ICC_p = mets["ICC_p"],
           CCC = mets["CCC"], CCC_lcl = mets["CCC_lcl"], CCC_ucl = mets["CCC_ucl"],
           d_peak_age = dpk, n_points = length(aligned$age))
  })
}

res_all <- do_one_split("overall")
if (compute_sex_strata) {
  res_all <- bind_rows(res_all, do_one_split("Female"), do_one_split("Male"))
}

# ----------------------------
# 4. SAVE OUTPUT TABLE
# ----------------------------
csv_path <- file.path(out_dir, "Pointwise_Metrics_Pearson_MAPE_ICC_CCC_Dpeak.csv")
write_csv(res_all, csv_path)
message("✅ Saved metrics table: ", csv_path)

# ----------------------------
# 5. VISUALIZATION
# ----------------------------
res_long <- res_all %>%
  filter(split == "overall") %>%
  select(Hemisphere, ROI, pearson_r, MAPE_refB, ICC, CCC, d_peak_age) %>%
  pivot_longer(cols = c(pearson_r, MAPE_refB, ICC, CCC, d_peak_age),
               names_to = "metric", values_to = "value")

metric_labels <- c(
  pearson_r = "Pearson r",
  MAPE_refB = "MAPE (ref = B)",
  ICC       = "ICC (2,1)",
  CCC       = "CCC (Lin)",
  d_peak_age = "ΔPeak age (years)"
)

p_box <- ggplot(res_long, aes(x = metric, y = value)) +
  geom_boxplot(outlier.shape = NA, fill = "#90CAF9", alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  scale_x_discrete(labels = metric_labels) +
  labs(x = NULL, y = NULL,
       title = paste0("Pointwise consistency: ", method_A_name, " vs ", method_B_name),
       subtitle = "Across hippocampal subfields (overall)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

pdf_path <- file.path(out_dir, "Pointwise_Metrics_Boxplot.pdf")
ggsave(pdf_path, p_box, width = 8, height = 5, device = cairo_pdf)
message("✅ Saved figure: ", pdf_path)

# ----------------------------
# 6. SUMMARY STATISTICS
# ----------------------------
mean_dpeak <- mean(res_all$d_peak_age, na.rm = TRUE)
median_dpeak <- median(res_all$d_peak_age, na.rm = TRUE)
range_dpeak <- range(res_all$d_peak_age, na.rm = TRUE)
cat("Mean Δpeak age:", mean_dpeak, "\n")
cat("Median Δpeak age:", median_dpeak, "\n")
cat("Range Δpeak age:", range_dpeak, "\n")

message("🎯 All analyses completed successfully.")

# ============================================================
# Multivariate SVR (RBF) with 10-fold CV
# Baseline & Longitudinal cognitive prediction in one script
# - Baseline: same-visit features -> same-visit outcome
# - Longitudinal: baseline features + follow-up time -> follow-up outcome
# - Evaluation: Pearson r (overall + per baseline Dx) with FDR
# GitHub-friendly paths
# ============================================================

# ==== Clean & Packages ====
rm(list = ls()); gc()
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr/ggplot2/tidyr/readr
  library(readxl)
  library(caret)
})

set.seed(2025_0812)

# ==== Paths ====
input_path <- "data/merged_from_tbl1.xlsx"     # put your xlsx here
out_root   <- "outputs/svr_cognitive"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ==== Config ====
id_col    <- "baseline_for_follow_up"          # subject linking column
time_var  <- "follow_up_time_days"             # >0 => follow-up
diag_keep <- c("MCI","AD","CSVD")              # baseline groups kept for longitudinal plots
target    <- "MMSE"                            # change to "MOCA" if needed

hippo_subfields <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)
base_covars <- c("Age","Sex","Diagnosis")     # as baseline covariates

# ==== Helpers ====
numize <- function(x) suppressWarnings(as.numeric(as.character(x)))
norm_sex <- function(x){
  z <- tolower(as.character(x))
  dplyr::case_when(
    z %in% c("male","m","1","男")   ~ 1,
    z %in% c("female","f","0","女") ~ 0,
    TRUE ~ NA_real_
  )
}
norm_diag <- function(x){
  y <- as.character(x)
  y[y %in% c("SVD","svd")] <- "CSVD"
  factor(y)
}
impute_med <- function(v){
  if (!is.numeric(v)) return(v)
  if (all(!is.finite(v))) return(rep(0, length(v)))
  v[!is.finite(v)] <- median(v[is.finite(v)], na.rm = TRUE)
  v
}

# ==== Read ====
df0 <- read_excel(input_path)

# Standardize key cols if present
for (v in c("Age","MOCA","MMSE", time_var)) if (v %in% names(df0)) df0[[v]] <- numize(df0[[v]])
if ("Sex" %in% names(df0))       df0$Sex       <- norm_sex(df0$Sex)
if ("Diagnosis" %in% names(df0)) df0$Diagnosis <- norm_diag(df0$Diagnosis)

# Optional: days -> years
# df0[[time_var]] <- df0[[time_var]] / 365.25

# ==== Identify baseline rows ====
is_baseline <- rep(FALSE, nrow(df0))
if ("baseline" %in% names(df0))
  is_baseline <- is_baseline | (tolower(as.character(df0$baseline)) == "baseline")
if (time_var %in% names(df0))
  is_baseline <- is_baseline | (!is.na(df0[[time_var]]) & df0[[time_var]] == 0)

# ------------------------------------------------------------------
# A) BASELINE MODEL  (same-visit features -> same-visit outcome)
# ------------------------------------------------------------------
run_baseline <- function(df) {
  out_dir <- file.path(out_root, "baseline"); dir.create(out_dir, showWarnings = FALSE)
  plot_dir <- file.path(out_dir, "plots");    dir.create(plot_dir, showWarnings = FALSE)
  
  # 选择一次就诊的数据（若同一受试者多行，可自行定义：这里简单取每个 ID 的最早一行）
  if (id_col %in% names(df)) {
    df_b <- df %>%
      filter(!is.na(.data[[id_col]]), .data[[id_col]] != "") %>%
      arrange(.data[[id_col]], dplyr::across(dplyr::any_of(time_var))) %>%
      group_by(.data[[id_col]]) %>% slice(1) %>% ungroup()
  } else {
    df_b <- df
  }
  
  # 建模表：同次就诊的 Age/Sex/Diagnosis + hippocampal + target
  cols_need <- c(base_covars, hippo_subfields, target)
  dat <- df_b %>%
    select(dplyr::anyOf(cols_need)) %>%
    filter(is.finite(.data[[target]]))
  
  # 缺失填补
  num_feats <- setdiff(intersect(cols_need, names(dat)), "Diagnosis")
  dat[num_feats] <- lapply(dat[num_feats], impute_med)
  
  # 拿 Diagnosis 作分组（不直接作为数值特征进入模型；如需编码，可 one-hot 后加入）
  group_vec <- as.character(dat$Diagnosis)
  x_cols <- setdiff(num_feats, target)
  
  # 10-fold CV SVR (RBF)
  cv_ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")
  tune_grid <- expand.grid(C = c(0.1,1,5,10),
                           sigma = c(0.01,0.05,0.1,0.2))
  svr_fit <- train(
    x = dat[, x_cols, drop = FALSE],
    y = dat[[target]],
    method = "svmRadial",
    trControl = cv_ctrl,
    preProcess = c("center","scale"),
    tuneGrid = tune_grid,
    metric   = "RMSE"
  )
  
  # OOF
  oof <- svr_fit$pred %>%
    filter(C == svr_fit$bestTune$C, sigma == svr_fit$bestTune$sigma) %>%
    arrange(rowIndex) %>%
    transmute(RowID = rowIndex, Predicted = pred)
  
  res <- tibble(RowID = seq_len(nrow(dat)),
                Actual = dat[[target]],
                Group  = factor(group_vec)) %>%
    left_join(oof, by = "RowID")
  
  # 评估
  cor_all <- suppressWarnings(cor.test(res$Actual, res$Predicted, method = "pearson"))
  groups  <- intersect(levels(res$Group), unique(as.character(res$Group)))
  raw_p <- sapply(groups, function(g){
    sub <- dplyr::filter(res, Group == g)
    if (nrow(sub) < 3) return(NA_real_)
    suppressWarnings(cor.test(sub$Actual, sub$Predicted, method = "pearson")$p.value)
  })
  names(raw_p) <- groups
  pFDRs <- p.adjust(raw_p, method = "fdr")
  
  # 保存
  readr::write_csv(res, file.path(out_dir, paste0("OOF_", target, "_baseline.csv")))
  saveRDS(svr_fit, file.path(out_dir, paste0("svr_fit_", target, "_baseline.rds")))
  
  # 作图
  color_map <- c("MCI"="#F39B7F", "AD"="#00A087", "CSVD"="#BC80BD", "PD"="#E64B35")
  for (g in groups) {
    sub <- dplyr::filter(res, Group == g)
    if (nrow(sub) < 3) next
    ct <- suppressWarnings(cor.test(sub$Actual, sub$Predicted, method = "pearson", conf.level = 0.95))
    pFDR_val <- pFDRs[g]
    p_txt <- if (is.na(pFDR_val)) "pFDR = NA"
    else if (pFDR_val < 0.001) "pFDR < 0.001"
    else sprintf("pFDR = %.3f", pFDR_val)
    ttl <- sprintf("Baseline %s\nr = %.2f [%.2f, %.2f], %s",
                   g, unname(ct$estimate), ct$conf.int[1], ct$conf.int[2], p_txt)
    
    p <- ggplot(sub, aes(Actual, Predicted)) +
      geom_point(color = color_map[g %||% "MCI"], alpha = .85, size = 3.6) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, color = "red", fill = "grey85") +
      labs(title = ttl, x = paste("Observed", target), y = paste("Predicted", target)) +
      coord_fixed(1) + theme_bw(base_size = 13) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 14, face = "bold", hjust = .5))
    ggsave(file.path(plot_dir, sprintf("SVR_%s_baseline_%s.pdf", target, g)),
           p, width = 5.5, height = 5.5, device = grDevices::cairo_pdf)
  }
  
  cat(sprintf("\n[BASELINE] %s: Overall Pearson r (OOF) = %.3f (p=%.3g); N=%d\n",
              target, unname(cor_all$estimate), cor_all$p.value, nrow(res)))
}

# ------------------------------------------------------------------
# B) LONGITUDINAL MODEL  (baseline feats + time -> follow-up outcome)
# ------------------------------------------------------------------
run_longitudinal <- function(df) {
  out_dir <- file.path(out_root, "longitudinal"); dir.create(out_dir, showWarnings = FALSE)
  plot_dir <- file.path(out_dir, "plots");        dir.create(plot_dir, showWarnings = FALSE)
  
  stopifnot(id_col %in% names(df))
  
  # baseline表
  is_baseline <- rep(FALSE, nrow(df))
  if ("baseline" %in% names(df))
    is_baseline <- is_baseline | (tolower(as.character(df$baseline)) == "baseline")
  if (time_var %in% names(df))
    is_baseline <- is_baseline | (!is.na(df[[time_var]]) & df[[time_var]] == 0)
  
  baseline_df <- df %>%
    filter(is_baseline, !is.na(.data[[id_col]]), .data[[id_col]] != "") %>%
    arrange(.data[[id_col]], dplyr::across(dplyr::any_of(time_var))) %>%
    group_by(.data[[id_col]]) %>% slice(1) %>% ungroup()
  
  baseline_feats <- baseline_df %>%
    select(dplyr::any_of(c(id_col, base_covars, hippo_subfields))) %>%
    mutate(Age = numize(Age), Sex = norm_sex(Sex), Diagnosis = norm_diag(Diagnosis)) %>%
    rename_with(~ paste0(.x, "_BL"), dplyr::any_of(c(base_covars, hippo_subfields)))
  
  # follow-up表
  follow_df <- df %>%
    filter(!is.na(.data[[time_var]]), .data[[time_var]] > 0,
           !is.na(.data[[id_col]]), .data[[id_col]] != "") %>%
    mutate(MOCA = numize(MOCA), MMSE = numize(MMSE)) %>%
    select(dplyr::any_of(c(id_col, time_var, "MOCA","MMSE")))
  
  dat_long <- follow_df %>%
    left_join(baseline_feats, by = id_col) %>%
    filter(is.na(Diagnosis_BL) | as.character(Diagnosis_BL) %in% diag_keep)
  
  feat_cols <- c(paste0(base_covars, "_BL"), paste0(hippo_subfields, "_BL"), time_var)
  dat <- dat_long %>%
    transmute(
      y = .data[[target]],
      Diagnosis_BL = if ("Diagnosis_BL" %in% names(dat_long)) factor(dat_long$Diagnosis_BL) else factor(NA),
      dplyr::across(dplyr::all_of(feat_cols))
    ) %>%
    filter(is.finite(y))
  
  num_feats <- intersect(feat_cols, names(dat))
  dat[num_feats] <- lapply(dat[num_feats], impute_med)
  dat <- dat %>% filter(complete.cases(dplyr::across(dplyr::all_of(num_feats))))
  
  # 10-fold CV SVR
  cv_ctrl <- trainControl(method = "cv", number = 10, savePredictions = "final")
  tune_grid <- expand.grid(C = c(0.1,1,5,10),
                           sigma = c(0.01,0.05,0.1,0.2))
  
  x_cols <- setdiff(num_feats, c("Diagnosis_BL"))
  svr_fit <- train(
    x = dat[, x_cols, drop = FALSE],
    y = dat$y,
    method = "svmRadial",
    trControl = cv_ctrl,
    preProcess = c("center","scale"),
    tuneGrid  = tune_grid,
    metric    = "RMSE"
  )
  
  oof <- svr_fit$pred %>%
    filter(C == svr_fit$bestTune$C, sigma == svr_fit$bestTune$sigma) %>%
    arrange(rowIndex) %>%
    transmute(RowID = rowIndex, Predicted = pred)
  
  res <- dat %>%
    mutate(RowID = dplyr::row_number()) %>%
    select(RowID, Actual = y, Group = Diagnosis_BL) %>%
    left_join(oof, by = "RowID")
  
  # 评估
  cor_all <- suppressWarnings(cor.test(res$Actual, res$Predicted, method = "pearson"))
  groups  <- intersect(diag_keep, unique(as.character(res$Group)))
  raw_p <- sapply(groups, function(g){
    sub <- dplyr::filter(res, Group == g)
    if (nrow(sub) < 3) return(NA_real_)
    suppressWarnings(cor.test(sub$Actual, sub$Predicted, method = "pearson")$p.value)
  })
  names(raw_p) <- groups
  pFDRs <- p.adjust(raw_p, method = "fdr")
  
  # 保存
  readr::write_csv(res, file.path(out_dir, paste0("OOF_", target, "_longitudinal.csv")))
  saveRDS(svr_fit, file.path(out_dir, paste0("svr_fit_", target, "_longitudinal.rds")))
  
  # 作图
  color_map <- c("MCI"="#F39B7F", "AD"="#00A087", "CSVD"="#BC80BD", "PD"="#E64B35")
  for (g in groups) {
    sub <- dplyr::filter(res, Group == g)
    if (nrow(sub) < 3) next
    ct <- suppressWarnings(cor.test(sub$Actual, sub$Predicted, method = "pearson", conf.level = 0.95))
    pFDR_val <- pFDRs[g]
    p_txt <- if (is.na(pFDR_val)) "pFDR = NA"
    else if (pFDR_val < 0.001) "pFDR < 0.001"
    else sprintf("pFDR = %.3f", pFDR_val)
    ttl <- sprintf("Longitudinal baseline=%s\nr = %.2f [%.2f, %.2f], %s",
                   g, unname(ct$estimate), ct$conf.int[1], ct$conf.int[2], p_txt)
    
    p <- ggplot(sub, aes(Actual, Predicted)) +
      geom_point(color = color_map[g %||% "MCI"], alpha = .85, size = 3.6) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, color = "red", fill = "grey85") +
      labs(title = ttl, x = paste("Observed", target), y = paste("Predicted", target)) +
      coord_fixed(1) + theme_bw(base_size = 13) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 14, face = "bold", hjust = .5))
    ggsave(file.path(plot_dir, sprintf("SVR_%s_longitudinal_%s.pdf", target, g)),
           p, width = 5.5, height = 5.5, device = grDevices::cairo_pdf)
  }
  
  cat(sprintf("\n[LONGITUDINAL] %s: Overall Pearson r (OOF) = %.3f (p=%.3g); N=%d\n",
              target, unname(cor_all$estimate), cor_all$p.value, nrow(res)))
}

# ==== Run both ====
run_baseline(df0)
run_longitudinal(df0)

cat("\nAll outputs written to: ", normalizePath(out_root, winslash = "/"), "\n", sep = "")

# ====================== Sex-stratified Welch + Cohen's d (HC - DIS, bootstrapped Welch ANOVA idea) ======================
rm(list = ls())

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(effectsize); library(rlang)
})

# -------- Config --------
input_path  <- "data.xlsx"
output_base <- "results"

restrict_age_by_disease <- TRUE   # 是否把 HC 限制在当前疾病的年龄范围
min_n_per_group <- 5              # 每组最小样本阈值（避免不稳）
B_boot <- 500                     # bootstrap 重采样次数（与方法描述一致）

# -------- Helpers --------
mean_ci_str <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2 || sd(x) == 0 || !is.finite(sd(x))) return(NA_character_)
  m <- mean(x); se <- sd(x) / sqrt(n); tcrit <- qt(0.975, df = n - 1)
  lo <- m - tcrit * se; hi <- m + tcrit * se
  sprintf(paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"), m, lo, hi)
}

# Cohen's d（方向：HC − DIS）
cohen_d_safe <- function(x_dis, y_hc) {
  df2 <- data.frame(
    value = c(y_hc, x_dis),                                         # 先 HC，后 DIS
    grp   = factor(c(rep("HC",  length(y_hc)), rep("DIS", length(x_dis))),
                   levels = c("HC","DIS"))                          # 关键：HC 在前
  )
  out <- tryCatch(
    effectsize::cohens_d(value ~ grp, data = df2, ci = 0.95, hedges.correction = FALSE),
    error = function(e) NULL
  )
  if (is.null(out)) return(c(NA_real_, NA_real_, NA_real_, NA_character_))
  d  <- out$Cohens_d[1]; lo <- out$CI_low[1]; hi <- out$CI_high[1]
  dstr <- if (all(is.finite(c(d,lo,hi)))) sprintf("%.3f [%.3f–%.3f]", d, lo, hi) else NA_character_
  c(d, lo, hi, dstr)
}

# ---- Bootstrapped non-parametric generalization of Welch’s one-way ANOVA (2-group case) ----
# 思路：在零假设下（两组均值相等），对每组做“中心化”（减去各自均值），
#      组内有放回重采样，保持组样本量与方差结构；每次计算 Welch t 统计量，
#      用 |t_boot| 的经验分布评估 |t_obs| 的双侧 p 值。
welch_boot_p <- function(x, y, B = 500, seed = NULL) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  if (!is.finite(sd(x)) || !is.finite(sd(y)) || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  if (!is.null(seed)) set.seed(seed)
  
  # 观测 Welch t
  t_obs <- tryCatch(as.numeric(t.test(x, y, var.equal = FALSE)$statistic),
                    error = function(e) NA_real_)
  if (!is.finite(t_obs)) return(NA_real_)
  
  # 组内中心化后的残差（零假设：均值相等）
  rx <- x - mean(x); ry <- y - mean(y)
  
  t_boot <- numeric(B)
  nx <- length(x); ny <- length(y)
  
  for (b in seq_len(B)) {
    xb <- sample(rx, size = nx, replace = TRUE)
    yb <- sample(ry, size = ny, replace = TRUE)
    # 在 H0 下两组均值都为 0，因此直接比较中心化后的重采样
    tb <- tryCatch(as.numeric(t.test(xb, yb, var.equal = FALSE)$statistic),
                   error = function(e) NA_real_)
    t_boot[b] <- tb
  }
  t_boot <- t_boot[is.finite(t_boot)]
  if (!length(t_boot)) return(NA_real_)
  
  # 双侧经验 p 值（加1修正避免0）
  (sum(abs(t_boot) >= abs(t_obs)) + 1) / (length(t_boot) + 1)
}

# -------- Data --------
df <- read_excel(input_path) %>%
  mutate(Diagnosis = recode(Diagnosis,
                            "MCI"="MCI","AD"="AD","PD"="PD","SVD"="CSVD","MS"="MS","HC"="HC"))

# 标准化 Sex -> Male / Female
if (is.numeric(df$Sex)) {
  df$Sex_std <- ifelse(df$Sex == 1, "Male", "Female")
} else {
  sx <- tolower(trimws(as.character(df$Sex)))
  df$Sex_std <- dplyr::case_when(
    sx %in% c("m","male","man","男","1") ~ "Male",
    sx %in% c("f","female","woman","女","0") ~ "Female",
    TRUE ~ NA_character_
  )
}
df$Sex_std <- factor(df$Sex_std, levels = c("Male","Female"))
df <- df %>% filter(!is.na(Sex_std)) %>% droplevels()

# -------- Settings --------
disease <- c("MCI","AD","PD","CSVD","MS")
subfields <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)

run_one_sex <- function(df_sex, sex_label, out_base) {
  results <- list()
  
  for (feature in subfields) for (dis in disease) {
    # 该性别内 HC + 当前疾病
    tmp <- df_sex %>%
      filter(Diagnosis %in% c("HC", dis)) %>%
      select(Age, Diagnosis, value = all_of(feature)) %>%
      mutate(value = as.numeric(value)) %>%
      filter(is.finite(Age) & is.finite(value)) %>%
      droplevels()
    
    if (!any(tmp$Diagnosis == "HC") || !any(tmp$Diagnosis == dis)) next
    
    # 可选：限制 HC 在疾病年龄范围内
    if (restrict_age_by_disease) {
      age_min <- min(tmp$Age[tmp$Diagnosis == dis], na.rm = TRUE)
      age_max <- max(tmp$Age[tmp$Diagnosis == dis], na.rm = TRUE)
      tmp <- tmp %>% filter(Age >= age_min, Age <= age_max)
    }
    
    # 划分两组
    x_dis <- tmp$value[tmp$Diagnosis == dis]
    y_hc  <- tmp$value[tmp$Diagnosis == "HC"]
    x_dis <- x_dis[is.finite(x_dis)]; y_hc <- y_hc[is.finite(y_hc)]
    
    n_dis <- length(x_dis); n_hc <- length(y_hc)
    if (min(n_dis, n_hc) < min_n_per_group) next
    
    mean_ci_DIS <- mean_ci_str(x_dis)
    mean_ci_HC  <- mean_ci_str(y_hc)
    
    # 变异性检查
    var_ok <- is.finite(sd(x_dis)) && is.finite(sd(y_hc)) && sd(x_dis) > 0 && sd(y_hc) > 0
    
    # ---- 关键替换：bootstrapped non-parametric generalization of Welch’s one-way ANOVA ----
    p_val <- if (var_ok) welch_boot_p(x_dis, y_hc, B = B_boot) else NA_real_
    
    # Cohen's d（方向：HC − DIS）
    d_parts <- if (var_ok) cohen_d_safe(x_dis, y_hc) else c(NA_real_, NA_real_, NA_real_, NA_character_)
    
    results[[length(results)+1]] <- data.frame(
      feature          = feature,
      disease          = dis,
      n_HC             = n_hc,
      n_dis            = n_dis,
      mean_ci_HC       = mean_ci_HC,
      mean_ci_dis      = mean_ci_DIS,
      p_value          = p_val,              # 这里就是 bootstrap Welch 的非参数 p 值
      cohens_d         = as.numeric(d_parts[1]),
      cohens_d_CI_low  = as.numeric(d_parts[2]),
      cohens_d_CI_high = as.numeric(d_parts[3]),
      cohens_d_str     = d_parts[4],
      stringsAsFactors = FALSE
    )
  }
  
  res <- dplyr::bind_rows(results)
  if (!nrow(res)) {
    warning("No results produced for ", sex_label, " (insufficient data).")
    return(invisible(NULL))
  }
  
  # 疾病内 BH-FDR
  res <- res %>%
    group_by(disease) %>%
    mutate(p_fdr = if (all(is.na(p_value))) NA_real_ else p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(disease, feature) %>%
    select(feature, disease, n_HC, n_dis, mean_ci_HC, mean_ci_dis,
           p_value, cohens_d, cohens_d_CI_low, cohens_d_CI_high, cohens_d_str, p_fdr)
  
  out_file <- paste0(out_base, "_", tolower(sex_label), ".csv")
  write.csv(res, out_file, row.names = FALSE)
  message("✅ Saved: ", out_file)
}

# -------- Run per sex --------
df_male   <- df %>% filter(Sex_std == "Male")
df_female <- df %>% filter(Sex_std == "Female")

run_one_sex(df_male,   "Male",   output_base)
run_one_sex(df_female, "Female", output_base)

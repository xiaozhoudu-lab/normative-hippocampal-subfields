# Associations between centile scores and cognitive measures:
# partial Spearman correlations controlling for age and sex, with BH-FDR.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(purrr);  library(writexl); library(ppcor)
})

# ---- 路径（相对路径，便于 GitHub 复现）----
input_path  <- "data/all.xlsx"            # 放示例或你自己的数据
output_dir  <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 读取 ----
df <- read_excel(input_path)

# ---- 变量 ----
# 注意：这里假定 subfields 列已经是“centile/quantile 或 deviation Z 分数”
subfields <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)
target <- "BVMT"   # 可换成 MoCA / MMSE / ……

# ---- 性别标准化为因子（Male/Female）----
if (!is.factor(df$Sex)) {
  sx <- tolower(trimws(as.character(df$Sex)))
  df$Sex <- dplyr::case_when(
    sx %in% c("m","male","man","男","1") ~ "Male",
    sx %in% c("f","female","woman","女","0") ~ "Female",
    TRUE ~ NA_character_
  ) |> factor(levels = c("Female","Male"))
}

# ---- 函数：偏 Spearman （对秩进行 Pearson 偏相关）----
partial_spearman <- function(x, y, age, sex) {
  # 秩变换 -> 在秩上做 Pearson 偏相关 ≈ partial Spearman
  ok <- is.finite(x) & is.finite(y) & is.finite(age) & !is.na(sex)
  if (sum(ok) < 10) return(c(NA_real_, NA_real_))  # 样本太少
  xr <- rank(x[ok]); yr <- rank(y[ok]); ar <- rank(age[ok])
  # sex 作为二元/多元哑变量处理：这里简单转 0/1（二分类）
  sx <- droplevels(sex[ok])
  if (nlevels(sx) == 2) {
    sx01 <- as.numeric(sx) - 1
    res <- ppcor::pcor.test(xr, yr, cbind(ar, sx01), method = "pearson")
  } else {
    # 多水平性别（少见），只控制 age
    res <- ppcor::pcor.test(xr, yr, ar, method = "pearson")
  }
  c(res$estimate, res$p.value)
}

# ---- 按诊断分组做“偏 Spearman + FDR”----
diagnosis_groups <- unique(df$Diagnosis)
all_results <- list()

for (diag in diagnosis_groups) {
  dfd <- df |> filter(Diagnosis == diag)
  
  # 只保留必要列，确保数值类型
  dfd <- dfd |>
    mutate(
      Age = as.numeric(Age),
      !!target := as.numeric(.data[[target]])
    )
  
  # 针对每个 subfield 做偏 Spearman
  cor_tbl <- map_dfr(subfields, function(region) {
    xr <- suppressWarnings(as.numeric(dfd[[region]]))
    out <- partial_spearman(
      x   = xr,
      y   = dfd[[target]],
      age = dfd$Age,
      sex = dfd$Sex
    )
    tibble(
      Diagnosis = diag,
      Region    = region,
      r_partial = out[1],
      p_value   = out[2]
    )
  })
  
  # FDR（组内多重校正）
  cor_tbl <- cor_tbl |>
    mutate(p_fdr = ifelse(all(is.na(p_value)), NA_real_, p.adjust(p_value, method = "BH")))
  
  all_results[[diag]] <- cor_tbl
}

final_results <- bind_rows(all_results)

# ---- 导出 ----
outfile <- file.path(output_dir, paste0("partial_spearman_", target, "_by_diagnosis.xlsx"))
writexl::write_xlsx(final_results, path = outfile)
message("Saved: ", outfile)

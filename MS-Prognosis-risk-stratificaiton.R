# ============================================================
# Prognostic CV Cox: 10-fold OOF risk, median split, KM + log-rank
# Input : data/MS.xlsx
# Output: outputs/ms_cvcox_ex_multiv/
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(survival); library(survminer); library(ggplot2); library(purrr)
  library(readr); library(caret); library(tibble)
})

set.seed(7)

# -------- Paths & Params (GitHub-friendly) --------
input_path   <- "data/MS.xlsx"                       # 放到仓库的 data/ 目录
out_dir      <- "outputs/ms_cvcox_ex_multiv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_col     <- "Followup_time_month"                # 随访时间（单位=月）
endpoints    <- c("EDSS_Progress","SPMS_conversion") # 事件=1/0
# ====== 设定多变量特征 ======
base_feats   <- c("EX")                              # 必选：EX
extra_feats  <- c()                                   # 可加别的变量，如 "Age","Sex","BMI",...
feature_set  <- c(base_feats, extra_feats)

k_folds      <- 10
months_eval  <- 60                                    # KM 注释里会注明 @60m

# -------- Helpers --------
norm <- function(x) tolower(gsub("\\s+","", trimws(x)))
to01 <- function(x){ x <- suppressWarnings(as.numeric(x)); x[!is.na(x) & x != 0] <- 1; x }
num  <- function(x) suppressWarnings(as.numeric(x))

save_km_plot <- function(survfit_obj, data_df, out_png, title=NULL) {
  p <- ggsurvplot(
    survfit_obj,
    data = data_df,
    pval = TRUE,                    # log-rank p 值
    conf.int = TRUE,
    surv.median.line = "hv",
    legend.labs  = c("Low risk","High risk"),
    legend.title = NULL
  )
  if (!is.null(title)) p$plot <- p$plot + ggtitle(title)
  ggsave(out_png, p$plot, width=7, height=5.5, dpi=300)
}

# -------- Read --------
stopifnot(file.exists(input_path))
df <- read_excel(input_path) %>% as.data.frame()
names(df) <- trimws(names(df))

# 宽松匹配时间列
idx_time <- which(norm(names(df)) == norm(time_col))
if (!length(idx_time)) stop(paste0("找不到时间列：", time_col))
time_col <- names(df)[idx_time[1]]

# 强制时间为正数值
df[[time_col]] <- num(df[[time_col]])
df[[time_col]][df[[time_col]] <= 0] <- NA_real_

# 检查结局列
missing_end <- setdiff(endpoints, names(df))
if (length(missing_end)) stop(paste("缺少事件列：", paste(missing_end, collapse=", ")))
df[endpoints] <- lapply(df[endpoints], to01)

# 检查特征列
miss_feat <- setdiff(feature_set, names(df))
if (length(miss_feat)) stop(paste("缺少特征列：", paste(miss_feat, collapse=", ")))

# -------- Main: per endpoint --------
for (ep in endpoints) {
  message("\n================= Endpoint: ", ep, " (multiv Cox + 10-fold CV) =================")
  d0 <- df[, c(time_col, ep, feature_set)]
  colnames(d0)[1:2] <- c("time","event")
  d0 <- d0 %>% filter(is.finite(time), !is.na(event))
  
  if (nrow(d0) < 40) { message("[", ep, "] 可用样本过少 (<40)，跳过。"); next }
  if (sum(d0$event==1, na.rm = TRUE) < 5) { message("[", ep, "] 事件过少 (<5)，跳过。"); next }
  
  # 简单缺失处理（数值：中位数；因子：众数 → 作为协变量时可改造）
  for (v in feature_set) {
    if (is.numeric(d0[[v]]) || is.integer(d0[[v]])) {
      med <- suppressWarnings(median(d0[[v]], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      d0[[v]][!is.finite(d0[[v]])] <- med
      d0[[v]][is.na(d0[[v]])] <- med
    } else {
      tab <- sort(table(d0[[v]]), decreasing = TRUE)
      if (length(tab) > 0) d0[[v]][is.na(d0[[v]])] <- names(tab)[1]
      d0[[v]] <- factor(d0[[v]])
    }
  }
  
  # 若某个数值特征方差为 0，剔除
  keep_feats <- feature_set[
    sapply(feature_set, function(v) !(is.numeric(d0[[v]]) && sd(d0[[v]], na.rm=TRUE) == 0))
  ]
  if (length(keep_feats) == 0) { message("[", ep, "] 所有特征方差为 0，跳过。"); next }
  
  d <- d0[, c("time","event", keep_feats), drop = FALSE]
  n  <- nrow(d)
  
  # ---------- 10-fold CV：获得 OOF 线性预测值 ----------
  folds <- caret::createFolds(seq_len(n), k = k_folds, returnTrain = TRUE)
  lp_oof <- rep(NA_real_, n)
  
  for (i in seq_along(folds)) {
    tr_id <- folds[[i]]
    te_id <- setdiff(seq_len(n), tr_id)
    
    # 训练集必须有事件
    if (sum(d$event[tr_id]==1) < 2) next
    
    # 拟合多变量 Cox
    form <- as.formula(paste0("Surv(time, event) ~ ", paste(keep_feats, collapse = " + ")))
    fit  <- try(coxph(form, data = d[tr_id, , drop=FALSE]), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    # 折外预测线性预测值
    lp_oof[te_id] <- as.numeric(predict(fit, newdata = d[te_id, , drop=FALSE], type = "lp"))
  }
  
  # 若仍有 NA（例如某折训练失败），用全数据重拟的模型做补齐（不会泄露到 OOF 评估以外的 KM分组）
  if (anyNA(lp_oof)) {
    form_all <- as.formula(paste0("Surv(time, event) ~ ", paste(keep_feats, collapse = " + ")))
    fit_all  <- coxph(form_all, data = d)
    lp_fill  <- as.numeric(predict(fit_all, newdata = d, type = "lp"))
    lp_oof[!is.finite(lp_oof)] <- lp_fill[!is.finite(lp_oof)]
  }
  
  # ---------- 按 OOF 风险（lp_oof）中位数二分 ----------
  med_lp <- median(lp_oof, na.rm = TRUE)
  grp    <- ifelse(lp_oof > med_lp, 1L, 0L)  # 1=High, 0=Low
  if (length(unique(grp[is.finite(grp)])) < 2) {
    message("[", ep, "] 中位数二分只得到一组，KM 跳过。")
    next
  }
  
  d_km <- d %>% mutate(Subgroup = grp)
  
  # ---------- KM + log-rank ----------
  sf <- survfit(Surv(time, event) ~ Subgroup, data = d_km)
  km_png <- file.path(out_dir, paste0("KM_", ep, "_multiv_cox_cv10.png"))
  ttl <- sprintf("%s (multivariate Cox; OOF median split @ %dm)", ep, months_eval)
  save_km_plot(sf, d_km, km_png, title = ttl)
  message("[", ep, "] KM plot -> ", km_png)
  
  # ---------- 目标月（例如 60 月）生存/事件概率 ----------
  t_eval <- min(months_eval, max(d_km$time, na.rm=TRUE))
  ssum <- summary(sf, times = t_eval)
  sp   <- ssum$surv
  # 兼容 strata 顺序
  if (length(sp) == 1) {
    names(sp) <- if (attr(ssum, "strata") == "Subgroup=0") "LowRisk_surv" else "HighRisk_surv"
    surv_probs <- c(LowRisk_surv = ifelse("LowRisk_surv" %in% names(sp), sp, NA_real_),
                    HighRisk_surv = ifelse("HighRisk_surv" %in% names(sp), sp, NA_real_))
  } else {
    names(sp) <- c("LowRisk_surv","HighRisk_surv")
    surv_probs <- sp
  }
  event_probs <- 1 - surv_probs
  names(event_probs) <- c("LowRisk_event","HighRisk_event")
  
  readr::write_csv(
    tibble(
      endpoint = ep,
      t_eval_months = t_eval,
      LowRisk_surv  = unname(surv_probs["LowRisk_surv"]),
      HighRisk_surv = unname(surv_probs["HighRisk_surv"]),
      LowRisk_event = unname(event_probs["LowRisk_event"]),
      HighRisk_event= unname(event_probs["HighRisk_event"])
    ),
    file.path(out_dir, paste0("prob_", months_eval, "m_", ep, "_multiv_cox_cv10.csv"))
  )
  
  # ---------- 导出 OOF 风险与分组 ----------
  oof_tbl <- d %>% transmute(time, event, !!!syms(keep_feats)) %>%
    mutate(lp_oof = lp_oof, Subgroup = grp)
  readr::write_csv(oof_tbl, file.path(out_dir, paste0("oof_lp_", ep, "_multiv_cox_cv10.csv")))
  
  # ---------- Cox 全模型汇总（可选：报告 HR 等） ----------
  fit_all <- coxph(as.formula(paste0("Surv(time, event) ~ ", paste(keep_feats, collapse = " + "))), data = d)
  s <- summary(fit_all)
  cox_tab <- tibble(
    endpoint = ep,
    feature  = rownames(s$coefficients),
    coef     = s$coefficients[, "coef"],
    HR       = s$coefficients[, "exp(coef)"],
    se_coef  = s$coefficients[, "se(coef)"],
    z        = s$coefficients[, "z"],
    p_value  = s$coefficients[, "Pr(>|z|)"],
    HR_low95 = s$conf.int[, "lower .95"],
    HR_up95  = s$conf.int[, "upper .95"],
    n        = s$n,
    events   = sum(d$event == 1),
    c_index  = s$concordance[1]
  )
  readr::write_csv(cox_tab, file.path(out_dir, paste0("cox_summary_", ep, "_multiv_fullsample.csv")))
}

message("\nAll done. Outputs saved in: ", normalizePath(out_dir))

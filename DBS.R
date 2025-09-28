# ============================================================
# SVM (RBF) with 10-fold CV producing OOF ROC
# Usage (optional args):
#   Rscript svm_oof_roc.R data/data.xlsx DBS outputs/
#   Rscript svm_oof_roc.R data/data.csv  DBS outputs/
# Defaults if args missing:
#   input_path = "data/data.xlsx", outcome_col = "DBS", out_dir = "outputs"
# ============================================================

# ---------- Packages ----------
need <- c("readxl","readr","dplyr","tidyr","stringr",
          "caret","kernlab","pROC","ggplot2")
inst <- setdiff(need, rownames(installed.packages()))
if (length(inst) > 0) install.packages(inst, dependencies = TRUE, quiet = TRUE)

suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr); library(stringr)
  library(caret);  library(kernlab); library(pROC); library(ggplot2)
})

# ---------- CLI args or defaults ----------
args <- commandArgs(trailingOnly = TRUE)
input_path  <- ifelse(length(args) >= 1, args[[1]], "data/data.xlsx")
outcome_col <- ifelse(length(args) >= 2, args[[2]], "DBS")
out_dir     <- ifelse(length(args) >= 3, args[[3]], "outputs")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# ---------- IO helpers ----------
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx","xls")) {
    readxl::read_excel(path)
  } else if (ext %in% c("csv","txt")) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

# ---------- Load data ----------
df <- read_any(input_path) %>% as.data.frame()
if (!outcome_col %in% names(df)) {
  stop("Outcome column '", outcome_col, "' not found. Available: ",
       paste(names(df), collapse=", "))
}

# ---------- Minimal cleaning ----------
# Outcome must be 0/1 numeric; we convert and check
df[[outcome_col]] <- suppressWarnings(as.numeric(df[[outcome_col]]))
if (!all(na.omit(df[[outcome_col]]) %in% c(0,1)) || any(is.na(df[[outcome_col]]))) {
  stop("'", outcome_col, "' must be strictly 0/1 with no NA.")
}

# ID (optional, for exporting OOF table)
if (!("ID" %in% names(df))) df$ID <- sprintf("row_%d", seq_len(nrow(df)))

# Select predictors = all except outcome
predictors <- setdiff(names(df), c(outcome_col))
# 常见无关列可剔除（按需）：时间戳、自由文本等
drop_guess <- c()  # c("Notes","Comment")
predictors <- setdiff(predictors, drop_guess)

dat <- df[, c("ID", predictors, outcome_col), drop = FALSE]

# 简单缺失值处理：数值→中位数；非数值→众数
impute_simple <- function(v){
  if (is.numeric(v)) {
    v[is.na(v)] <- median(v, na.rm = TRUE)
    v
  } else {
    v <- as.character(v)
    if (any(is.na(v))) {
      tab <- sort(table(v), decreasing = TRUE)
      fill <- if (length(tab)) names(tab)[1] else "UNK"
      v[is.na(v)] <- fill
    }
    factor(v)
  }
}
for (nm in predictors) dat[[nm]] <- impute_simple(dat[[nm]])

# caret 需要分类标签因子：定义正类为“Responder”=1
dat[[outcome_col]] <- factor(ifelse(dat[[outcome_col]] == 1, "Responder", "NonResponder"),
                             levels = c("NonResponder","Responder"))

# ---------- Outer 10-fold (OOF) ----------
folds <- createFolds(dat[[outcome_col]], k = 10, returnTrain = TRUE)

# Inner CV for tuning (5-fold); standardize in training pipeline
ctrl_inner <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE, summaryFunction = twoClassSummary,
  # 失衡可开启下采样（按需）：
  sampling = "down",
  savePredictions = "final"
)

# 存放折外概率（对“Responder”的概率）
oof_prob <- rep(NA_real_, nrow(dat))
names(oof_prob) <- dat$ID

# 仅特征矩阵（去掉 ID & outcome）
X_cols <- setdiff(colnames(dat), c("ID", outcome_col))

for (tr in folds) {
  te <- setdiff(seq_len(nrow(dat)), tr)
  
  # 训练本折模型（SVM-RBF，自动调参；指标 AUC）
  mdl <- caret::train(
    x = dat[tr, X_cols, drop = FALSE],
    y = dat[[outcome_col]][tr],
    method = "svmRadial",
    metric = "ROC",
    preProcess = c("center","scale"),
    tuneLength = 10,
    trControl = ctrl_inner
  )
  
  # 预测测试折概率（正类 Responder）
  pr <- predict(mdl, newdata = dat[te, X_cols, drop = FALSE], type = "prob")
  oof_prob[ dat$ID[te] ] <- pr$Responder
}

# ---------- ROC on OOF ----------
valid <- is.finite(oof_prob)
if (!any(valid)) stop("No valid OOF probabilities produced. Check data and features.")
roc_obj <- pROC::roc(response = dat[[outcome_col]][valid],
                     predictor = oof_prob[valid],
                     levels = c("NonResponder","Responder"),
                     ci = TRUE)

auc_val <- as.numeric(pROC::auc(roc_obj))
ci_val  <- as.numeric(roc_obj$ci)

# 保存 OOF 表
oof_tbl <- data.frame(
  ID = dat$ID,
  true_label = dat[[outcome_col]],
  prob_responder = oof_prob,
  stringsAsFactors = FALSE
)
write.csv(oof_tbl, file.path(out_dir, "OOF_predictions.csv"), row.names = FALSE)

# 保存 AUC 文本
cat(sprintf("SVM (RBF) 10-fold OOF AUC = %.3f [%.3f–%.3f]\n",
            auc_val, ci_val[1], ci_val[3]),
    file = file.path(out_dir, "AUC_OOF.txt"))

# 绘制并保存 ROC 图（PNG）
png(file.path(out_dir, "ROC_OOF.png"), width = 900, height = 700, res = 120)
plot(roc_obj,
     main = sprintf("OOF ROC (SVM-RBF, 10-fold)  AUC=%.3f [%.3f–%.3f]",
                    auc_val, ci_val[1], ci_val[3]),
     col = "#2C7FB8", lwd = 3)
abline(0,1,lty=2,col="grey60")
dev.off()

# 也可导出 ggplot 风格
gg <- ggroc(roc_obj, colour = "#2C7FB8", size = 1.2) +
  geom_abline(slope=1, intercept=1, linetype="dashed", colour="grey60") +
  theme_classic(base_size = 14) +
  labs(title = sprintf("OOF ROC (SVM-RBF, 10-fold)  AUC=%.3f [%.3f–%.3f]",
                       auc_val, ci_val[1], ci_val[3]),
       x = "1 - Specificity", y = "Sensitivity")
ggplot2::ggsave(file.path(out_dir, "ROC_OOF_ggplot.png"), gg, width = 7.5, height = 6, dpi = 150)

cat("✔ Done. Outputs saved in: ", normalizePath(out_dir, winslash = "/"), "\n", sep = "")

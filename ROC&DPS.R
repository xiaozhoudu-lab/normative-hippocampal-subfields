#!/usr/bin/env Rscript
# dps_svm_vsHC_matched.R
# Train disease-vs-HC SVM (RBF) models with 10-fold CV.
# In each outer fold, perform 1:1 nearest-neighbor matching on Age + Sex (training fold only).
# Produce out-of-fold (OOF) AUC and OOF DPS for all subjects.
# Outputs: AUC_vsHC_OOF.csv, DPS_OOF_vsHC_AllSubjects.csv, Fig_ROC_vsHC.pdf, DPS_Radar_AllModels_OOF.pdf

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl)
  library(caret); library(kernlab); library(pROC)
  library(ggplot2); library(grDevices); library(MatchIt)
})

# ------------------------- CLI args / Defaults -------------------------
args <- commandArgs(trailingOnly = TRUE)
# Usage examples:
# Rscript scripts/dps_svm_vsHC_matched.R data/all.xlsx outputs
# Rscript scripts/dps_svm_vsHC_matched.R data/all.csv   outputs

INPUT_FILE <- if (length(args) >= 1) args[1] else "data/all.xlsx"
OUTPUT_DIR <- if (length(args) >= 2) args[2] else "outputs"
RADAR_DIR  <- file.path(OUTPUT_DIR, "DPS_Radars")
SEED <- 123

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RADAR_DIR,  recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

# Colors for plots
col_vsHC <- c("MCI"="#F8766D","AD"="#00BFC4","PD"="#C77CFF","CSVD"="#7CAE00","MS"="#FF7F00")
model_colors <- col_vsHC
open_pdf <- function(file, w=8, h=8) grDevices::cairo_pdf(file, width=w, height=h, family="Helvetica")

# ------------------------- IO helpers -------------------------
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx","xls")) {
    readxl::read_excel(path)
  } else if (ext %in% c("csv","txt")) {
    read.csv(path, stringsAsFactors = FALSE)
  } else stop("Unsupported file extension: ", ext)
}

# ------------------------- Load data -------------------------
df <- read_any(INPUT_FILE) %>%
  mutate(
    Diagnosis = dplyr::recode(
      Diagnosis,
      "HC" = "HC", "MCI" = "MCI", "AD" = "AD", "PD" = "PD",
      "SVD" = "CSVD", "MS" = "MS", "AQP4Pos_NMOSD" = "NMOSD"
    )
  ) %>%
  filter(Diagnosis != "NMOSD")

# Create ID if missing
if (".SubjectID" %in% names(df)) {
  df$ID <- as.character(df$.SubjectID)
} else if (!("ID" %in% names(df))) {
  df$ID <- sprintf("row_%d", seq_len(nrow(df)))
} else {
  df$ID <- as.character(df$ID)
}
df$ID[is.na(df$ID) | df$ID == ""] <- sprintf("row_%d", which(is.na(df$ID) | df$ID == ""))
df$ID <- make.unique(df$ID, sep = "_")

# Normalize Sex to numeric (0/1) for matching
df <- df %>%
  mutate(
    Sex = case_when(
      Sex %in% c("Male","M",1,"1") ~ 1,
      Sex %in% c("Female","F",0,"0") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(Sex))
    )
  )

# Hippocampal subfields (24)
subfields <- c(
  "L.subiculum","L.CA1","L.presubiculum","L.molecular_layer_HP",
  "L.GC.ML.DG","L.CA3","L.CA4","L.Whole_hippocampus",
  "L.Hippocampal_tail","L.parasubiculum","L.fimbria","L.HATA",
  "R.subiculum","R.CA1","R.presubiculum","R.molecular_layer_HP",
  "R.GC.ML.DG","R.CA3","R.CA4","R.Whole_hippocampus",
  "R.fimbria","R.parasubiculum","R.Hippocampal_tail","R.HATA"
)

# Which diseases to model (vs HC)
diseases <- c("MCI","AD","PD","CSVD","MS")

# Keep needed cols; coerce features to numeric
need_cols <- unique(c("ID","Age","Sex","Diagnosis", subfields))
df <- df[, intersect(need_cols, names(df))]
df <- df %>% mutate(across(any_of(subfields), ~ suppressWarnings(as.numeric(.x))))

# ------------------------- Outer 10-fold CV on IDs -------------------------
folds <- createFolds(df$ID, k = 10, returnTrain = TRUE)

# Inner CV for model tuning
ctrl_inner <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE, summaryFunction = twoClassSummary,
  sampling = "down"  # class imbalance: down-sample majority in training resamples
)

# ------------------------- Core routine (Age+Sex matching in TRAIN fold) -------------------------
fit_one_vsHC_get_OOF <- function(pos_label, feats) {
  oof_task_prob <- rep(NA_real_, nrow(df))  # OOF prob for AUC (disease vs HC only)
  oof_all_prob  <- rep(NA_real_, nrow(df))  # OOF prob for all subjects (DPS)
  
  for (tr in folds) {
    te <- setdiff(seq_len(nrow(df)), tr)
    train_df <- df[tr, ]
    test_df  <- df[te, ]
    
    # Training subset: pos_label vs HC
    train_bin0 <- train_df %>%
      filter(Diagnosis %in% c(pos_label, "HC")) %>%
      select(all_of(c("ID","Diagnosis","Age","Sex", feats))) %>%
      filter(is.finite(Age), is.finite(Sex))
    
    # If either class missing, skip
    if (!any(train_bin0$Diagnosis == pos_label) || !any(train_bin0$Diagnosis == "HC")) next
    
    # Limit HC in training to disease age range
    age_min <- min(train_bin0$Age[train_bin0$Diagnosis == pos_label], na.rm = TRUE)
    age_max <- max(train_bin0$Age[train_bin0$Diagnosis == pos_label], na.rm = TRUE)
    train_bin0 <- train_bin0 %>%
      filter(!(Diagnosis == "HC" & (Age < age_min | Age > age_max)))
    
    # Age + Sex nearest-neighbor 1:1 matching (TRAIN only)
    # Make disease the first level so twoClassSummary treats it as "event"
    train_bin0$Diagnosis <- factor(train_bin0$Diagnosis, levels = c(pos_label, "HC"))
    m.out <- try(
      MatchIt::matchit(Diagnosis ~ Age + Sex, data = train_bin0,
                       method = "nearest", ratio = 1),
      silent = TRUE
    )
    if (inherits(m.out, "try-error")) next
    train_matched <- MatchIt::match.data(m.out)
    
    # Remove near-zero variance features & rows with missing features
    nzv <- nearZeroVar(train_matched[, feats, drop = FALSE])
    feats_use <- if (length(nzv)) feats[-nzv] else feats
    if (!length(feats_use)) next
    
    train_matched <- train_matched %>%
      filter(if_all(all_of(feats_use), ~ is.finite(.))) %>%
      droplevels()
    
    if (nrow(train_matched) < 10 || nlevels(train_matched$Diagnosis) < 2) next
    
    # Train SVM (RBF) with inner CV on matched training set
    mdl <- caret::train(
      x = train_matched[, feats_use, drop = FALSE],
      y = train_matched$Diagnosis,  # levels: c(pos_label, "HC")
      method = "svmRadial",
      preProcess = c("center","scale"),
      metric = "ROC",
      tuneLength = 10,
      trControl = ctrl_inner
    )
    
    # (1) OOF for AUC: test fold's disease + HC only
    test_task <- test_df %>%
      filter(Diagnosis %in% c(pos_label, "HC")) %>%
      select(all_of(c("ID","Diagnosis", feats)))
    if (nrow(test_task)) {
      ok <- stats::complete.cases(test_task[, feats_use, drop = FALSE])
      if (any(ok)) {
        pr <- predict(mdl, newdata = test_task[ok, feats_use, drop = FALSE], type = "prob")
        idx <- match(test_task$ID[ok], df$ID)
        # Probability for the positive class = pos_label
        oof_task_prob[idx] <- pr[[pos_label]]
      }
    }
    
    # (2) OOF for DPS: test fold's ALL subjects
    ok2 <- stats::complete.cases(test_df[, feats_use, drop = FALSE])
    if (any(ok2)) {
      pr2 <- predict(mdl, newdata = test_df[ok2, feats_use, drop = FALSE], type = "prob")
      idx2 <- match(test_df$ID[ok2], df$ID)
      oof_all_prob[idx2] <- pr2[[pos_label]]
    }
  }
  
  list(prob_task = oof_task_prob, prob_all = oof_all_prob)
}

# ------------------------- Run per disease -------------------------
oof_list <- lapply(diseases, fit_one_vsHC_get_OOF, feats = subfields)
names(oof_list) <- diseases

# ------------------------- AUC table (OOF, disease vs HC) -------------------------
auc_df <- do.call(rbind, lapply(diseases, function(dx) {
  idx <- which(df$Diagnosis %in% c(dx, "HC") & is.finite(oof_list[[dx]]$prob_task))
  if (!length(idx)) return(data.frame(Disease = dx, AUC = NA, CI_L = NA, CI_U = NA))
  roc_obj <- pROC::roc(response = factor(df$Diagnosis[idx], levels = c("HC", dx)),
                       predictor = oof_list[[dx]]$prob_task[idx], ci = TRUE)
  data.frame(
    Disease = dx,
    AUC  = as.numeric(pROC::auc(roc_obj)),
    CI_L = as.numeric(pROC::ci.auc(roc_obj)[1]),
    CI_U = as.numeric(pROC::ci.auc(roc_obj)[3])
  )
}))
write.csv(auc_df, file.path(OUTPUT_DIR, "AUC_vsHC_OOF.csv"), row.names = FALSE)

# ------------------------- DPS wide table (OOF for ALL subjects) -------------------------
dps_oof <- data.frame(ID = df$ID, Diagnosis = df$Diagnosis)
for (dx in diseases) dps_oof[[dx]] <- oof_list[[dx]]$prob_all
write.csv(dps_oof, file.path(OUTPUT_DIR, "DPS_OOF_vsHC_AllSubjects.csv"), row.names = FALSE)

# ------------------------- ROC overlay PDF -------------------------
roc_list <- list()
for (dx in diseases) {
  idx <- which(df$Diagnosis %in% c(dx, "HC") & is.finite(oof_list[[dx]]$prob_task))
  if (length(idx)) {
    roc_list[[dx]] <- pROC::roc(
      response  = factor(df$Diagnosis[idx], levels = c("HC", dx)),
      predictor = oof_list[[dx]]$prob_task[idx], ci = TRUE
    )
  } else roc_list[[dx]] <- NULL
}

roc_pdf <- file.path(OUTPUT_DIR, "Fig_ROC_vsHC.pdf")
open_pdf(roc_pdf, 8, 8); on.exit(grDevices::dev.off(), add = TRUE)
op <- par(mar = c(5, 6.5, 1, 2) + 0.1, xaxs = "i", yaxs = "i", family = "Helvetica")

valid_dx <- names(roc_list)[vapply(roc_list, Negate(is.null), logical(1))]
if (length(valid_dx) == 0) {
  plot.new(); text(0.5,0.5,"No valid ROC to plot")
} else {
  dx1  <- valid_dx[1]; roc1 <- roc_list[[dx1]]
  fpr1 <- 1 - roc1$specificities; tpr1 <- roc1$sensitivities; o1 <- order(fpr1)
  plot(fpr1[o1], tpr1[o1], type="l", lwd=3, col=model_colors[dx1],
       xlim=c(0,1), ylim=c(0,1), asp=1,
       xlab="False Positive Rate (1 - Specificity)",
       ylab="True Positive Rate (Sensitivity)",
       cex.lab=1.8, cex.axis=1.6, font.lab=2, font.axis=2)
  abline(0,1,lty=2)
  if (length(valid_dx) > 1) {
    for (dx in valid_dx[-1]) {
      roci <- roc_list[[dx]]
      fpr  <- 1 - roci$specificities; tpr <- roci$sensitivities; o <- order(fpr)
      lines(fpr[o], tpr[o], lwd=3, col=model_colors[dx])
    }
  }
  legend_text <- vapply(valid_dx, function(dx){
    aucv <- as.numeric(pROC::auc(roc_list[[dx]]))
    civ  <- as.numeric(pROC::ci.auc(roc_list[[dx]]))
    sprintf("%s  AUC=%.3f [%.3f–%.3f]", dx, aucv, civ[1], civ[3])
  }, character(1))
  legend("bottomright", legend = legend_text, col = model_colors[valid_dx],
         lwd = 3, cex = 1.4, text.font = 2, bty = "n")
}
par(op); grDevices::dev.off()

# ------------------------- DPS radar (median OOF per group) -------------------------
diag_levels <- c("HC","MCI","AD","PD","CSVD","MS")
model_names <- diseases
radar_source <- dps_oof

long <- radar_source %>%
  dplyr::select(Diagnosis, dplyr::any_of(model_names)) %>%
  tidyr::pivot_longer(-Diagnosis, names_to = "Model", values_to = "Value")

agg <- long %>%
  dplyr::group_by(Model, Diagnosis) %>%
  dplyr::summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::complete(Model = model_names, Diagnosis = diag_levels, fill = list(Value = 0)) %>%
  dplyr::mutate(
    Diagnosis = factor(Diagnosis, levels = diag_levels),
    Model     = factor(Model,     levels = model_names),
    Value     = pmin(pmax(Value, 0), 1)
  ) %>%
  dplyr::arrange(Model, Diagnosis)

K <- length(diag_levels)
theta_vec <- (pi/2) - 2*pi*(0:(K-1))/K
theta_map <- setNames(as.numeric(theta_vec), diag_levels)

vals_xy <- agg %>%
  dplyr::mutate(
    theta = theta_map[as.character(Diagnosis)],
    x = Value * cos(theta),
    y = Value * sin(theta)
  )

vals_closed <- vals_xy %>%
  dplyr::group_by(Model) %>%
  dplyr::arrange(Diagnosis, .by_group = TRUE) %>%
  dplyr::do(dplyr::bind_rows(., dplyr::slice(., 1))) %>%
  dplyr::ungroup()

labs_df <- data.frame(
  Diagnosis = factor(diag_levels, levels = diag_levels),
  x = 1.10 * cos(theta_vec),
  y = 1.10 * sin(theta_vec),
  lab = diag_levels
)
circle_df <- function(r, n=720){ t <- seq(0, 2*pi, length.out = n); data.frame(x=r*cos(t), y=r*sin(t)) }

p_radar <- ggplot() +
  geom_path(data = circle_df(1.00), aes(x, y), linetype="dashed", color="grey75", linewidth=0.7) +
  geom_path(data = circle_df(0.50), aes(x, y), linetype="dashed", color="grey85", linewidth=0.6) +
  geom_path(data = circle_df(0.10), aes(x, y), linetype="dashed", color="grey92", linewidth=0.5) +
  geom_segment(aes(x=0, y=0, xend=0, yend=1), color="grey80", linewidth=0.6) +
  geom_path(data = vals_closed, aes(x=x, y=y, group=Model, color=Model), linewidth=1.4, lineend="round") +
  geom_point(data = vals_xy, aes(x=x, y=y, color=Model), shape=21, fill="white", stroke=1.1, size=2.8) +
  geom_text(data = labs_df, aes(x=x, y=y, label=lab), color="grey30", size=4, fontface="bold") +
  annotate("text", x=0, y=0.02, label="0%",   size=3.6, color="grey50") +
  annotate("text", x=0, y=0.50, label="50%",  size=3.6, color="grey50", vjust=-0.3) +
  annotate("text", x=0, y=1.00, label="100%", size=3.6, color="grey50", vjust=-0.3) +
  scale_color_manual(values = model_colors, name = "Models",
                     guide = guide_legend(override.aes = list(shape = 21, fill = "white", size = 3))) +
  coord_fixed() +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.background = element_rect(fill="white", color=NA),
    plot.background  = element_rect(fill="white", color=NA),
    panel.grid       = element_blank(),
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(face="bold"),
    legend.text      = element_text(face="bold"),
    plot.title       = element_text(size=16, face="bold", hjust=0.5),
    plot.margin      = margin(8,8,8,8)
  ) +
  ggtitle("DPS (OOF) from disease-vs-HC models")

pdf_path <- file.path(RADAR_DIR, "DPS_Radar_AllModels_OOF.pdf")
ggsave(pdf_path, p_radar, width=8, height=8, units="in",
       device=grDevices::cairo_pdf, bg="white")

message("✅ Outputs:\n- ", file.path(OUTPUT_DIR, "AUC_vsHC_OOF.csv"),
        "\n- ", file.path(OUTPUT_DIR, "DPS_OOF_vsHC_AllSubjects.csv"),
        "\n- ", file.path(OUTPUT_DIR, "Fig_ROC_vsHC.pdf"),
        "\n- ", file.path(RADAR_DIR,  "DPS_Radar_AllModels_OOF.pdf"))

# ============================================================
# 11_Sex_stratified_sensitivity.R
# Within-disease comparison of male and female hippocampal
# normative percentile scores (24 bilateral regions).
# ============================================================

rm(list = ls())

# Usage (recommended; output from 02_Apply_normative_model.R):
# Rscript 11_Sex_stratified_sensitivity.R \
#   <normative_scores_long.csv> <clinical_data.csv> <output_dir> \
#   [diagnoses] [minimum_per_sex] [diagnosis_column]
#
# The script also accepts the legacy wide .xlsx file used to generate the
# original sex_difference_outputs.  In that case, use "-" for clinical_data:
# Rscript 11_Sex_stratified_sensitivity.R \
#   <legacy_centile_workbook.xlsx> - <output_dir>
#
# diagnoses is a comma-separated list.  Defaults: MCI,AD,PD,CSVD,MS.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript 11_Sex_stratified_sensitivity.R ",
    "<normative_scores_long.csv|legacy_centile_workbook.xlsx> ",
    "<clinical_data.csv|-> <output_dir> [diagnoses] ",
    "[minimum_per_sex] [diagnosis_column]"
  )
}

score_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
clinical_arg <- args[[2]]
output_dir <- args[[3]]
disease_arg <- if (length(args) >= 4) args[[4]] else "MCI,AD,PD,CSVD,MS"
minimum_per_sex <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
diagnosis_column <- if (length(args) >= 6) args[[6]] else "Diagnosis"

if (!is.finite(minimum_per_sex) || minimum_per_sex < 2) {
  stop("minimum_per_sex must be an integer of at least 2")
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

normalize_sex <- function(x) {
  raw <- tolower(safe_chr(x))
  out <- rep(NA_character_, length(raw))
  out[raw %in% c("male", "m", "man", "1")] <- "Male"
  out[raw %in% c("female", "f", "woman", "0", "2")] <- "Female"
  out
}

normalize_diagnosis <- function(x) {
  out <- toupper(safe_chr(x))
  out[out %in% c("SVD", "CSVD")] <- "CSVD"
  out
}

make_match_id <- function(df) {
  if ("match_id" %in% names(df)) return(safe_chr(df$match_id))
  candidates <- list(
    c("Freesufer_Path2", "Freesufer_Path3"),
    c("Freesurfer_Path2", "Freesurfer_Path3")
  )
  for (pair in candidates) {
    if (all(pair %in% names(df))) {
      return(paste0(safe_chr(df[[pair[[1]]]]), safe_chr(df[[pair[[2]]]])))
    }
  }
  stop(
    "Data require match_id or matching FreeSurfer path columns ",
    "(Freesufer/Freesurfer_Path2 and Path3)."
  )
}

normalize_feature_name <- function(x) {
  out <- safe_chr(x)
  out <- sub("_centile$", "", out, ignore.case = TRUE)
  out <- sub("^L_", "left_", out)
  out <- sub("^R_", "right_", out)
  out <- sub("^lh[._]", "left_", out, ignore.case = TRUE)
  out <- sub("^rh[._]", "right_", out, ignore.case = TRUE)
  out <- gsub("GC_ML_DG", "GC.ML.DG", out, fixed = TRUE)
  out <- sub("_molecular_layer$", "_molecular_layer_HP", out)
  out
}

feature_label <- function(feature_name) {
  x <- normalize_feature_name(feature_name)
  hemi <- ifelse(startsWith(tolower(x), "left_"), "L.", "R.")
  roi <- sub("^(left|right)_", "", x, ignore.case = TRUE)
  roi <- gsub("Whole_hippocampus", "hippocampus", roi, fixed = TRUE)
  roi <- gsub("GC.ML.DG", "GC-DG", roi, fixed = TRUE)
  roi <- gsub("GC_ML_DG", "GC-DG", roi, fixed = TRUE)
  roi <- gsub("molecular_layer_HP", "ML", roi, fixed = TRUE)
  roi <- gsub("molecular_layer", "ML", roi, fixed = TRUE)
  roi <- gsub("Hippocampal_tail", "tail", roi, fixed = TRUE)
  paste0(hemi, roi)
}

feature_order <- function(feature_name) {
  x <- normalize_feature_name(feature_name)
  hemi <- ifelse(startsWith(tolower(x), "left_"), 0L, 1L)
  roi <- sub("^(left|right)_", "", x, ignore.case = TRUE)
  roi <- gsub("GC_ML_DG", "GC.ML.DG", roi, fixed = TRUE)
  roi <- gsub("molecular_layer$", "molecular_layer_HP", roi)
  order <- c(
    "Whole_hippocampus", "CA1", "CA3", "CA4", "GC.ML.DG",
    "molecular_layer_HP", "presubiculum", "subiculum",
    "parasubiculum", "HATA", "fimbria", "Hippocampal_tail"
  )
  index <- match(roi, order)
  index[is.na(index)] <- 99L
  hemi * 100L + index
}

load_legacy_workbook <- function(path, diagnosis_column) {
  selected <- NULL
  selected_sheet <- NULL
  for (sheet in excel_sheets(path)) {
    candidate <- try(as.data.frame(read_excel(path, sheet = sheet)), silent = TRUE)
    if (inherits(candidate, "try-error")) next
    centile_columns <- grep("_centile$", names(candidate), value = TRUE, ignore.case = TRUE)
    if (all(c(diagnosis_column, "Sex") %in% names(candidate)) &&
        length(centile_columns) > 0) {
      selected <- candidate
      selected_sheet <- sheet
      break
    }
  }
  if (is.null(selected)) {
    stop(
      "No worksheet contains Diagnosis, Sex, and columns ending in _centile: ",
      path
    )
  }

  selected$match_id <- if ("match_id" %in% names(selected)) {
    safe_chr(selected$match_id)
  } else {
    paste0("legacy_row_", seq_len(nrow(selected)))
  }

  centile_columns <- grep("_centile$", names(selected), value = TRUE, ignore.case = TRUE)
  long_parts <- lapply(centile_columns, function(column) {
    data.frame(
      match_id = selected$match_id,
      Diagnosis = selected[[diagnosis_column]],
      Sex = selected$Sex,
      feature_name = normalize_feature_name(column),
      Quant_score = safe_num(selected[[column]]),
      stringsAsFactors = FALSE
    )
  })

  message("Using legacy worksheet: ", selected_sheet)
  bind_rows(long_parts)
}

load_normative_scores <- function(score_file, clinical_arg, diagnosis_column) {
  extension <- tolower(tools::file_ext(score_file))
  if (extension %in% c("xlsx", "xls")) {
    return(load_legacy_workbook(score_file, diagnosis_column))
  }
  if (extension != "csv") stop("score_file must be .csv, .xlsx, or .xls")

  scores <- read.csv(score_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("match_id", "hemisphere", "ROI", "Quant_score")
  missing_required <- setdiff(required, names(scores))
  if (length(missing_required) > 0) {
    stop(
      "CSV input must be normative_scores_long.csv from ",
      "02_Apply_normative_model.R. Missing: ",
      paste(missing_required, collapse = ", ")
    )
  }
  if (clinical_arg == "-") {
    stop("clinical_data.csv is required with normative_scores_long.csv")
  }

  clinical_file <- normalizePath(clinical_arg, winslash = "/", mustWork = TRUE)
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!(diagnosis_column %in% names(clinical))) {
    stop("Clinical data do not contain diagnosis column: ", diagnosis_column)
  }
  clinical$match_id <- make_match_id(clinical)
  clinical <- clinical[clinical$match_id != "", , drop = FALSE]
  if (anyDuplicated(clinical$match_id)) {
    stop("Clinical data contain duplicated match_id values")
  }

  clinical_keep <- data.frame(
    match_id = clinical$match_id,
    Diagnosis_clinical = clinical[[diagnosis_column]],
    stringsAsFactors = FALSE
  )
  if ("Sex" %in% names(clinical)) clinical_keep$Sex_clinical <- clinical$Sex

  scores$match_id <- safe_chr(scores$match_id)
  scores$feature_name <- paste0(tolower(safe_chr(scores$hemisphere)), "_", safe_chr(scores$ROI))
  scores$Quant_score <- safe_num(scores$Quant_score)
  if (!("Sex" %in% names(scores))) scores$Sex <- NA_character_

  merged <- inner_join(scores, clinical_keep, by = "match_id")
  if (nrow(merged) == 0) stop("No score rows matched the clinical data")

  sex_scores <- normalize_sex(merged$Sex)
  if ("Sex_clinical" %in% names(merged)) {
    sex_clinical <- normalize_sex(merged$Sex_clinical)
    disagreement <- !is.na(sex_scores) & !is.na(sex_clinical) & sex_scores != sex_clinical
    if (any(disagreement)) {
      stop("Sex differs between normative scores and clinical data for matched participants")
    }
    sex_scores[is.na(sex_scores)] <- sex_clinical[is.na(sex_scores)]
  }

  data.frame(
    match_id = merged$match_id,
    Diagnosis = merged$Diagnosis_clinical,
    Sex = sex_scores,
    feature_name = normalize_feature_name(merged$feature_name),
    Quant_score = merged$Quant_score,
    stringsAsFactors = FALSE
  )
}

analysis_data <- load_normative_scores(score_file, clinical_arg, diagnosis_column)
analysis_data$Diagnosis <- normalize_diagnosis(analysis_data$Diagnosis)
analysis_data$Sex <- normalize_sex(analysis_data$Sex)
analysis_data$feature_name <- normalize_feature_name(analysis_data$feature_name)
analysis_data$Quant_score <- safe_num(analysis_data$Quant_score)

diseases <- normalize_diagnosis(strsplit(disease_arg, ",", fixed = TRUE)[[1]])
diseases <- unique(diseases[diseases != ""])
if (length(diseases) == 0) stop("No valid diagnoses were supplied")

analysis_data <- analysis_data[
  analysis_data$Diagnosis %in% diseases &
    analysis_data$Sex %in% c("Male", "Female") &
    is.finite(analysis_data$Quant_score),
  , drop = FALSE
]
if (nrow(analysis_data) == 0) stop("No eligible disease centile records were found")

outside_probability <- analysis_data$Quant_score < -1e-10 | analysis_data$Quant_score > 1 + 1e-10
if (any(outside_probability)) {
  stop("Quant_score contains values outside the expected percentile range [0, 1]")
}
analysis_data$Quant_score <- pmin(pmax(analysis_data$Quant_score, 0), 1)

duplicate_key <- duplicated(analysis_data[, c("match_id", "feature_name")])
if (any(duplicate_key)) {
  stop("Duplicated match_id + feature_name rows were found in the score input")
}

expected_features <- as.vector(outer(
  c("left", "right"),
  c(
    "Whole_hippocampus", "CA1", "CA3", "CA4", "GC.ML.DG",
    "molecular_layer_HP", "presubiculum", "subiculum",
    "parasubiculum", "HATA", "fimbria", "Hippocampal_tail"
  ),
  paste,
  sep = "_"
))
observed_features <- unique(analysis_data$feature_name)
missing_features <- setdiff(expected_features, observed_features)
if (length(missing_features) > 0) {
  warning(
    "The input does not contain all 24 expected bilateral regions. Missing: ",
    paste(missing_features, collapse = ", ")
  )
}

stats_rows <- list()
skip_rows <- list()
row_index <- 0L
skip_index <- 0L

for (disease in diseases) {
  disease_data <- analysis_data[analysis_data$Diagnosis == disease, , drop = FALSE]
  for (feature in sort(unique(disease_data$feature_name))) {
    feature_data <- disease_data[disease_data$feature_name == feature, , drop = FALSE]
    male <- feature_data$Quant_score[feature_data$Sex == "Male"]
    female <- feature_data$Quant_score[feature_data$Sex == "Female"]
    male <- male[is.finite(male)]
    female <- female[is.finite(female)]
    n_male <- length(male)
    n_female <- length(female)

    if (n_male < minimum_per_sex || n_female < minimum_per_sex) {
      skip_index <- skip_index + 1L
      skip_rows[[skip_index]] <- data.frame(
        Disease = disease,
        feature_name = feature,
        n_male = n_male,
        n_female = n_female,
        reason = "fewer than minimum_per_sex observations",
        stringsAsFactors = FALSE
      )
      next
    }

    male_sd <- stats::sd(male)
    female_sd <- stats::sd(female)
    pooled_variance <- (
      (n_male - 1) * male_sd^2 + (n_female - 1) * female_sd^2
    ) / (n_male + n_female - 2)
    cohen_d <- if (is.finite(pooled_variance) && pooled_variance > 0) {
      (mean(male) - mean(female)) / sqrt(pooled_variance)
    } else {
      NA_real_
    }
    se_d <- if (is.finite(cohen_d)) {
      sqrt(
        (n_male + n_female) / (n_male * n_female) +
          cohen_d^2 / (2 * (n_male + n_female - 2))
      )
    } else {
      NA_real_
    }
    welch <- try(stats::t.test(male, female, var.equal = FALSE), silent = TRUE)

    row_index <- row_index + 1L
    stats_rows[[row_index]] <- data.frame(
      Diagnosis_code = disease,
      Disease = disease,
      Subfield = feature_label(feature),
      Column = feature,
      Subfield_order = feature_order(feature),
      Hemisphere = ifelse(startsWith(feature, "left_"), "Left", "Right"),
      n_male = n_male,
      n_female = n_female,
      male_mean = mean(male),
      female_mean = mean(female),
      male_sd = male_sd,
      female_sd = female_sd,
      cohen_d_male_minus_female = cohen_d,
      d_ci_low = cohen_d - 1.96 * se_d,
      d_ci_high = cohen_d + 1.96 * se_d,
      welch_t = if (inherits(welch, "try-error")) NA_real_ else unname(welch$statistic),
      welch_df = if (inherits(welch, "try-error")) NA_real_ else unname(welch$parameter),
      p_value = if (inherits(welch, "try-error")) NA_real_ else welch$p.value,
      stringsAsFactors = FALSE
    )
  }
}

if (length(stats_rows) == 0) {
  stop("No disease-region comparison met minimum_per_sex")
}

statistics <- bind_rows(stats_rows)
statistics$q_fdr_within_disease <- NA_real_
for (disease in unique(statistics$Disease)) {
  index <- which(statistics$Disease == disease)
  statistics$q_fdr_within_disease[index] <- stats::p.adjust(
    statistics$p_value[index], method = "BH"
  )
}
statistics$q_fdr <- statistics$q_fdr_within_disease
statistics$abs_d <- abs(statistics$cohen_d_male_minus_female)
statistics$direction <- ifelse(
  statistics$cohen_d_male_minus_female >= 0,
  "Male higher", "Female higher"
)
statistics <- statistics[order(statistics$Disease, -statistics$abs_d), , drop = FALSE]

write.csv(
  statistics,
  file.path(output_dir, "sex_difference_stats_all_diseases_subfields.csv"),
  row.names = FALSE
)

if (length(skip_rows) > 0) {
  write.csv(
    bind_rows(skip_rows),
    file.path(output_dir, "sex_difference_skip_log.csv"),
    row.names = FALSE
  )
}

participant_data <- unique(analysis_data[, c("match_id", "Diagnosis", "Sex")])
count_long <- participant_data %>%
  count(Diagnosis, Sex, name = "n")
count_grid <- expand.grid(
  Disease = diseases,
  Sex = c("Female", "Male"),
  stringsAsFactors = FALSE
)
count_long <- left_join(
  count_grid,
  count_long,
  by = c("Disease" = "Diagnosis", "Sex" = "Sex")
)
count_long$n[is.na(count_long$n)] <- 0L
sample_counts <- data.frame(
  Disease = diseases,
  Female = count_long$n[match(paste(diseases, "Female"), paste(count_long$Disease, count_long$Sex))],
  Male = count_long$n[match(paste(diseases, "Male"), paste(count_long$Disease, count_long$Sex))],
  stringsAsFactors = FALSE
)
write.csv(
  sample_counts,
  file.path(output_dir, "sex_difference_sample_counts.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    setting = c(
      "score_file", "clinical_file", "diagnosis_column", "diagnoses",
      "minimum_per_sex", "effect_direction", "test", "multiplicity"
    ),
    value = c(
      score_file, clinical_arg, diagnosis_column, paste(diseases, collapse = ","),
      minimum_per_sex, "Male minus Female", "two-sided Welch t-test",
      "BH-FDR separately within each disease across regions"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "sex_difference_settings.csv"),
  row.names = FALSE
)

significance_symbol <- function(q) {
  ifelse(
    is.na(q), "",
    ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))
  )
}

disease_colors <- c(
  AD = "#76A9B7", CSVD = "#4A567A", MCI = "#B4504A",
  MS = "#CC8C78", PD = "#5C947A"
)

plot_data <- statistics
plot_data$Disease <- factor(plot_data$Disease, levels = diseases)
subfield_levels <- unique(
  plot_data$Subfield[order(plot_data$Subfield_order)]
)
plot_data$Subfield <- factor(plot_data$Subfield, levels = subfield_levels)
plot_data$tile_label <- paste0(
  sprintf("%.2f", plot_data$cohen_d_male_minus_female),
  significance_symbol(plot_data$q_fdr)
)

heatmap_plot <- ggplot(
  plot_data,
  aes(x = Subfield, y = Disease, fill = cohen_d_male_minus_female)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = tile_label), size = 2.7) +
  scale_fill_gradient2(
    low = "#31678B", mid = "#F6F6F6", high = "#B24641",
    midpoint = 0, name = "Cohen's d"
  ) +
  labs(
    title = "Sex differences in hippocampal normative percentiles",
    subtitle = "Male minus Female; stars indicate within-disease BH-FDR",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "sex_difference_all_subfields_heatmap.png"),
  heatmap_plot, width = 14, height = 5.5, dpi = 300
)
ggsave(
  file.path(output_dir, "sex_difference_all_subfields_heatmap_editable.pdf"),
  heatmap_plot, width = 14, height = 5.5, device = cairo_pdf
)

significant <- statistics[
  is.finite(statistics$q_fdr) & statistics$q_fdr < 0.05,
  , drop = FALSE
]

if (nrow(significant) > 0) {
  split_significant <- split(
    significant,
    interaction(significant$Disease, significant$Hemisphere, drop = TRUE)
  )
  strongest <- bind_rows(lapply(split_significant, function(x) {
    head(x[order(-x$abs_d), , drop = FALSE], 3)
  }))
  strongest$effect_label <- paste(strongest$Disease, strongest$Subfield, sep = " - ")
  strongest$effect_label <- factor(
    strongest$effect_label,
    levels = rev(unique(strongest$effect_label[order(strongest$Disease, strongest$Hemisphere)]))
  )

  forest_plot <- ggplot(
    strongest,
    aes(x = cohen_d_male_minus_female, y = effect_label, color = Disease)
  ) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.5) +
    geom_segment(
      aes(x = d_ci_low, xend = d_ci_high, yend = effect_label),
      linewidth = 0.8
    ) +
    geom_point(size = 2.4) +
    facet_wrap(~Hemisphere, scales = "free_y") +
    scale_color_manual(values = disease_colors, drop = FALSE) +
    labs(
      title = "Strongest FDR-significant sex differences",
      subtitle = "Up to three effects per disease and hemisphere",
      x = "Cohen's d (Male - Female)", y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  box_keys <- head(significant[order(-significant$abs_d), , drop = FALSE], 24)
  box_data <- inner_join(
    analysis_data,
    box_keys[, c("Diagnosis_code", "Column", "Subfield", "q_fdr")],
    by = c("Diagnosis" = "Diagnosis_code", "feature_name" = "Column")
  )
  box_data$panel_label <- paste0(
    box_data$Diagnosis, " - ", box_data$Subfield, " ",
    significance_symbol(box_data$q_fdr)
  )
  box_data$Sex <- factor(box_data$Sex, levels = c("Male", "Female"))
  box_plot <- ggplot(box_data, aes(x = Sex, y = Quant_score, fill = Sex)) +
    geom_boxplot(width = 0.65, outlier.alpha = 0.25) +
    scale_fill_manual(values = c(Male = "#76A9B7", Female = "#B4504A")) +
    facet_wrap(~panel_label, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "Centile distributions for significant disease-region effects",
      x = NULL, y = "Normative percentile"
    ) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(size = 7)
    )
} else {
  forest_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No FDR-significant sex differences") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() +
    labs(title = "Strongest FDR-significant sex differences")
  box_plot <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No significant effects to display") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void()
}

save_summary_figure <- function(path, type = c("png", "pdf")) {
  type <- match.arg(type)
  if (type == "png") {
    grDevices::png(path, width = 2400, height = 2600, res = 300)
  } else {
    grDevices::cairo_pdf(path, width = 10, height = 11)
  }
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2, ncol = 1,
    heights = grid::unit(c(0.43, 0.57), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(forest_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(box_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grDevices::dev.off()
}

save_summary_figure(
  file.path(output_dir, "sex_difference_summary_figure.png"), "png"
)
save_summary_figure(
  file.path(output_dir, "sex_difference_summary_figure_editable.pdf"), "pdf"
)

message(
  "Sex-stratified sensitivity analysis completed: ",
  nrow(statistics), " disease-region comparisons; ",
  sum(statistics$q_fdr < 0.05, na.rm = TRUE), " within-disease FDR-significant."
)

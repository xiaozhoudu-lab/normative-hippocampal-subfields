# Hippocampal Subfield Normative Modeling

This repository accompanies the study:

**Mapping hippocampal subfield heterogeneity across neurological disorders by normative modeling**

It contains:

1. the analysis scripts used for normative modeling and downstream analyses;
2. public normative-reference resources for bilateral whole-hippocampal and hippocampal-subfield volumes; and
3. a lightweight function for calculating centile scores in new data.

The workflow constructs lifespan normative references, calculates individualized Z scores and centile scores, evaluates disease-related deviations, and performs clinical, prognostic, and sensitivity analyses.

Centile scores are represented on a 0–1 scale. Some legacy output names use the term `percentile` for compatibility; in this repository, `centile score` and `percentile score` refer to the same quantity.

## Repository structure

```text
.
├── scripts/
│   ├── 01_Normative_model_fit.R
│   ├── 02_Apply_normative_model.R
│   ├── 03_New_site_calibration.R
│   ├── 04_Group_comparison.R
│   ├── 05_CMD.R
│   ├── 06_Classification_DPS.R
│   ├── 07_Cognitive_prediction.R
│   ├── 08_Progression_analysis.R
│   ├── 09_ICV_sensitivity.R
│   ├── 10_ComBatGAM_sensitivity.R
│   └── 11_Sex_stratified_sensitivity.R
├── R/
│   ├── 100.common-variables.r
│   ├── 101.common-functions.r
│   ├── 102.gamlss-recode.r
│   ├── 300.variables.r
│   ├── 301.functions.r
│   └── ZZZ_function.R
├── Check-and-install-packages.R
├── Hippo_Subfields_Normative_PublicRelease.zip
├── Hippo_Subfields_Normative_PublicRelease_Age4_85_ICV.zip
├── LICENSE
└── README.md
```

Scripts are numbered in their recommended execution order. Sensitivity scripts can be run after the primary normative models and individualized scores have been generated.

## Workflow

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_Normative_model_fit.R` | Fit the primary lifespan normative models and generate normative trajectories. |
| 2 | `02_Apply_normative_model.R` | Apply fitted models and calculate individualized Z and centile scores. |
| 3 | `03_New_site_calibration.R` | Calibrate pretrained models for previously unseen sites or scanners. |
| 4 | `04_Group_comparison.R` | Compare disease groups with 1:1 matched healthy controls. |
| 5 | `05_CMD.R` | Calculate centile Mahalanobis distance (CMD). |
| 6 | `06_Classification_DPS.R` | Perform disease classification and calculate disease propensity scores (DPS). |
| 7 | `07_Cognitive_prediction.R` | Perform cognitive association and prediction analyses. |
| 8 | `08_Progression_analysis.R` | Perform longitudinal progression and survival analyses. |
| 9 | `09_ICV_sensitivity.R` | Fit ICV-adjusted normative models. |
| 10 | `10_ComBatGAM_sensitivity.R` | Evaluate ComBat-GAM harmonization and empirical site-scale effects. |
| 11 | `11_Sex_stratified_sensitivity.R` | Evaluate sex-stratified disease-related centile patterns. |

## Public normative-reference resources

Two public release packages are provided:

- `Hippo_Subfields_Normative_PublicRelease.zip`: primary hippocampal and hippocampal-subfield normative reference;
- `Hippo_Subfields_Normative_PublicRelease_Age4_85_ICV.zip`: ICV-adjusted reference for ages 4–85 years.

Each archive contains a normative quantile grid, release metadata, a lightweight scoring function, example input data, instructions, and reference figures. Internal fitted model objects are not distributed because they may retain information linked to study sites or participant-level data.

After extracting an archive, new observations can be scored using:

```r
source("R/score_normative_centile.R")

score_normative_file(
  input_csv = "examples/example_subjects.csv",
  output_csv = "examples/example_subjects_scored.csv",
  grid_path = "data/hippo_subfields_normative_quantile_grid.csv.gz"
)
```

The input file must contain `hemisphere`, `roi`, `sex`, `age`, and `volume`. Scoring is based on interpolation within the released quantile grid.

## Script descriptions

### `01_Normative_model_fit.R`

Fits GAMLSS normative models for bilateral whole-hippocampal and hippocampal-subfield volumes. The models include nonlinear age effects, sex, and site/scanner effects and generate normative trajectories and participant-level Z and centile scores.

### `02_Apply_normative_model.R`

Applies the fitted normative models to a cohort or individual observations. Input data must use compatible volume units, FreeSurfer definitions, sex coding, and regional names.

### `03_New_site_calibration.R`

Estimates new-site location and scale offsets using eligible healthy controls from sites or scanners absent from model training. It then applies the calibrated model to participants from the corresponding site and reports calibration status for downstream interpretation.

### `04_Group_comparison.R`

Compares each disease group with 1:1 healthy controls matched exactly on site/scanner and sex and matched to the nearest available age without replacement. Inference uses a Welch statistic with labels permuted within matched pairs; two-sided empirical P values are based on 5,000 permutations, and regional P values are BH-FDR corrected within each disease group. CMD can be evaluated using the same matching framework as a separate test family.

### `05_CMD.R`

Calculates CMD directly from the 24-dimensional bilateral hippocampal centile profile. The healthy-control mean vector and covariance matrix define the reference space, and higher CMD indicates greater multivariate deviation from the healthy-control profile. The fitted reference can be saved and applied to an independent cohort.

### `06_Classification_DPS.R`

Randomly divides each HC-versus-disease dataset into 70% training and 30% independent held-out test partitions, followed by separate 1:1 age- and sex-matching within each partition. Regularized logistic regression and single-hidden-layer neural-network models are compared by 5-fold cross-validation in the training data; the selected model is evaluated once in the held-out test data. Held-out probabilities are used to calculate ROC AUC, DeLong 95% confidence intervals, and DPS.

The script requires one unique clinical record and one value per participant-feature combination. Duplicate participant identifiers or participant-feature rows trigger an error rather than being silently retained as independent observations.

### `07_Cognitive_prediction.R`

Includes partial Spearman associations with age and sex adjustment, cross-sectional cognitive prediction, and longitudinal cognitive prediction. Regional association P values are BH-FDR corrected within each analysis-group-by-cognitive-outcome family; CMD–cognition tests form a separate family within each analysis group.

Cognitive prediction uses radial-basis-function support vector regression with outer 10-fold cross-validation. Standardization and model tuning are performed using training-fold data only, and performance is quantified by Pearson correlation across pooled held-out predictions. Repeated observations from the same participant remain in the same fold.

### `08_Progression_analysis.R`

Performs leakage-controlled prognostic modeling for the HC/MCI and MS endpoints. Participants with incomplete survival outcomes or model predictors are excluded before cross-validation; no imputation is performed. LASSO Cox feature selection and fitting occur within the training folds, and the fitted coefficients generate risk scores only for the corresponding held-out folds.

Harrell's concordance index is calculated from pooled out-of-fold risk scores. Kaplan–Meier visualization and two-sided log-rank testing use high- and low-risk groups defined by the median pooled out-of-fold risk score.

### `09_ICV_sensitivity.R`

Repeats the normative-model workflow with intracranial volume included as a covariate in the GAMLSS location parameter. Individualized scores use each participant's ICV, whereas reference trajectories are evaluated at the healthy-control median ICV.

### `10_ComBatGAM_sensitivity.R`

Evaluates whether site/scanner variation materially affects normative estimates. The first analysis applies ComBat-GAM to healthy-control hippocampal volumes while preserving nonlinear age and sex effects before refitting the normative models. The second analysis perturbs the primary model's scale parameter using empirical healthy-control site-specific variability and compares the resulting centiles with the primary estimates.

### `11_Sex_stratified_sensitivity.R`

Compares male and female centile scores within each disease using two-sided Welch t-tests. Effect sizes are reported as Male-minus-Female Cohen's d, with BH-FDR correction performed separately within each disease across regional tests.

## Shared modeling functions

Files in `R/` contain the shared configuration, modified GAMLSS utilities, model-fitting functions, parameter extraction, prediction, and project-specific centile-scoring functions required by the numbered scripts. They are dependencies rather than standalone analyses.

## Expected input data

The scripts require appropriately formatted participant-level clinical information and FreeSurfer-derived hippocampal volume tables. Depending on the analysis, typical fields include:

- participant or scan identifier;
- age and sex;
- site or scanner identifier;
- diagnosis;
- bilateral whole-hippocampal and subfield volumes;
- intracranial volume for the ICV sensitivity analysis; and
- cognitive scores or longitudinal outcomes.

Example data, where supplied, demonstrate file structure and code execution only and must not be used for scientific inference.

## Software and required R packages

The primary analyses were performed in R 4.3.1. Hippocampal segmentation and extraction of MRI-derived measures used FreeSurfer 7.3.2. Detailed software and package information is provided in Supplementary Method S6.

Install the packages required by scripts 01–11 from the repository root:

```r
source("Check-and-install-packages.R")
```

The direct dependencies are `readxl`, `dplyr`, `stringr`, `reshape2`, `gamlss`, `gamlss.add`, `gamlss.dist`, `foreach`, `doParallel`, `caret`, `glmnet`, `nnet`, `pROC`, `kernlab`, `survival`, `mgcv`, `ComBatFamily`, and `ggplot2`. All are installed from CRAN except `ComBatFamily`, which is obtained from its [official GitHub repository](https://github.com/andy1764/ComBatFamily).

## Data availability and privacy

Original individual-level MRI and clinical datasets are not redistributed because access and redistribution are governed by institutional approvals, participant consent, and data-use agreements. Third-party cohort data remain available through their original access procedures. No personally identifiable information or original individual-level clinical data are included in this repository.

## Reproducibility

The repository provides the principal modeling and statistical procedures described in the manuscript. Because the original datasets cannot be redistributed, users must combine the scripts with appropriately approved data formatted according to the supplied examples. Cohort definitions, MRI processing, quality control, statistical design, external validation, and sensitivity analyses are described in the manuscript and Supplementary Methods.

## Citation

If you use these scripts or normative-reference resources, please cite:

**Mapping hippocampal subfield heterogeneity across neurological disorders by normative modeling.**  
Accepted for publication in *Nature Communications*.

The full citation and DOI will be added after online publication.

## License

This project is released under the MIT License. See `LICENSE` for details.

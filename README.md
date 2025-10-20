# Charting hippocampal subfield normative modeling and its clinical applications

## Table of contents

1. [Introduction](#Introduction)
2. [Datasets](#Datasets)
3. [Models](#Models)
4. [Scripts](#Scripts)
   4.1 [Required R packages and installation](#Check-and-install-packages.R)
   4.2 [Normative curve estimation and peak age determination](#Normative-model-fit.R)
   4.3 [Bootstrap analysis of normative curves](#Bootstrap-normative-mdoel-fit.R)
   4.4 [ICV adjustment and model comparison](#ICV.R)
   4.5 [Model calibration using new datasets](#Calibration-normative-model-using-new-dataset.R)
   4.6 [Applying normative models to disease datasets](#Disease-application-normative-model.R)
   4.7 [Applying normative models to individual-level data](#Individual-application-normative-model.R)
   4.8 [Statistical analyses of deviation scores across diseases](#Stastical-analysis-deviations-across-diseases.R)
   4.9 [Clinical downstream tasks](#Clinical-tasks)
5. [License](#License)

---

## 1. Introduction

This repository provides scripts for constructing **hippocampal subfield normative models** using *Generalized Additive Models for Location, Scale and Shape (GAMLSS)* and for applying these models to diverse clinical tasks, including disease classification, cognitive prediction, and prognosis stratification.

The repository includes:

* **Example datasets** for normative model fitting, calibration, and validation.
* **Trained normative models (RDS files)** that can be requested from the authors.
* **Scripts** for each stage of the analytical workflow.
* **Example outputs** for reference and validation.

---

## 2. Datasets

Example datasets are provided solely for demonstration and code testing purposes.
They include small synthetic samples that mimic the format of the original datasets.
These examples are **not suitable for scientific analyses**.

The full multi-center datasets used in published studies cannot be shared publicly due to data-use restrictions.
Researchers interested in collaboration or access to the trained normative models should contact the authors directly via email.

---

## 3. Models

The pre-trained **GAMLSS-based normative models** for hippocampal subfield volumes are available as `.rds` files.
These models can be applied to new subjects or cohorts to derive age-, sex-, and site-adjusted deviation scores (Z-scores and centiles).

> 🔹 **Note:**
> The complete RDS model files are not included in this repository due to data-sharing policies.
> Researchers can request access by contacting the corresponding author via email.

---

## 4. Scripts

The repository contains a complete, modular workflow for normative modeling and clinical applications.

### 4.1 Required R packages and installation

**`Check-and-install-packages.R`**
Automatically checks and installs all necessary R packages required for modeling, application, and visualization.

---

### 4.2 Normative curve estimation and peak age determination

**`Normative-model-fit.R`**
Fits normative lifespan models using GAMLSS based on hippocampal subfield volumes, estimates peak ages, and generates normative trajectories.

---

### 4.3 Bootstrap analysis of normative curves

**`Bootstrap-normative-mdoel-fit.R`**
Performs bootstrap resampling (default 1,000 iterations) to estimate confidence intervals for normative trajectories and peak ages.

---

### 4.4 ICV adjustment and model comparison

**`ICV.R`**
Implements intracranial volume (ICV)–adjusted normative models and evaluates consistency between adjusted and non-adjusted trajectories using the concordance correlation coefficient (CCC).

---

### 4.5 Model calibration using new datasets

**`Calibration-normative-model-using-new-dataset.R`**
Applies a recalibration procedure to adapt pre-trained normative models to new sites or scanners using control subjects.

---

### 4.6 Applying normative models to disease datasets

**`Disease-application-normative-model.R`**
Applies the trained normative models to multi-disease datasets to compute deviation scores (Z-scores and centiles) for each patient across all hippocampal subfields.

---

### 4.7 Applying normative models to individual-level data

**`Individual-application-normative-model.R`**
Computes deviation scores for a single individual or small cohort based on pre-trained normative models, enabling personalized brain health assessments.

---

### 4.8 Statistical analyses of deviation scores across diseases

**`Stastical-analysis-deviations-across-diseases.R`**
Performs Monte Carlo permutation tests for group-level comparisons between diseases and matched controls, estimating p-values and effect sizes (Cohen’s d) with FDR correction.

---

### 4.9 Clinical downstream tasks

**`ROC&DPS.R`**
Performs ROC analyses and computes disease propensity scores (DPS) using SVM-based classification of deviation profiles.

**`associations_centile_scores_cognition.R`**
Examines correlations between hippocampal centile scores and cognitive measures (partial Spearman correlation, FDR-corrected).

**`cognitive-score-prediciton.R`**
Uses multivariate SVR (RBF kernel) to predict cognitive performance from centile scores.

**`HC&MCI-risk-stratificaiton.R`**
Conducts longitudinal survival analysis to predict conversion from healthy control (HC) to mild cognitive impairment (MCI), and from MCI to Alzheimer’s disease (AD).

**`MS-Prognosis-risk-stratificaiton.R`**
Builds prognostic models for multiple sclerosis (MS) progression using LASSO–Cox regression and Kaplan–Meier survival analysis.

**`DBS.R`**
Implements an SVM classifier to identify responders vs. non-responders to deep brain stimulation (DBS) treatment in Parkinson’s disease.

---

## 5. License

**MIT License**

Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


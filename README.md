# normative-hippocampal-subfields
Charting hippocampal subfield normative modeling and its clinical applications

## Table of contents

1. [Introduction](#Introduction)
2. [Datasets](#Datasets)
   2.1. Example dataset for the normative reference estimation (#Dataset-norms)
   2.2. Example disease dataset for the normative reference application (#Dataset-diseases)
   2.3. Example individual data for the personalized application of normative reference (#Dataset-individual)
   2.4. Example new dataset for the site-calibration of normative reference (#Dataset-new)
   2.5. Supplementary/Extended datasets for a better interpretation of the normative reference application (#Dataset-extended)
3. [Normative models and DPS models](#Models)
   3.1 Normative models by GAMLSS based on hippocampal subfields (#GAMLSS)
   3.2 Disease-specific classification models for DPS estimation (#DPS)
4. [Scripts](#Scripts)
   4.1 Required R packages and installation (#Check-and-install-packages.R)
   4.2 Normative curve estimate and milestone determination (#Normative-model-fit.R)
   4.3 Bootstrap analysis of the estimated normative curves (#Bootstrap-normative-mdoel-fit.R)
   4.4 Apply normative model to disease datasets (#Disease-application-normative-model.R)
   4.5 Apply normative model to an individual case (#Individual-application-normative-model.R)
   4.6 Calibrate normative model by new dataset (#Calibration-normative-model-using-new-dataset.R)
   4.7 Statistical analyses of deviation score across diseases (#Stastical-analysis-deviations-across-diseases.R)
   4.8 Other clinical task related scripts (#ROC&DPS.R; #cognitive-score-prediciton.R; #DBS.R; #HC&MCI-risk-stratificaiton.R; #MS-Prognosis-risk-stratificaiton.R; #associations_centile_scores_cognition.R)
5. [License](#License)

---

## 1. Introduction

This repository contains codes for hippocampal subfield normative reference construction and their downstream clinical applications. The example datasets could be found in file **“Datasets”**; the normative models (references) could be found in file **“Models”**; the main scripts could be found in file **“Scripts”**; the test outputs by the main scripts using example datasets are provided in file **“Test_results”**.

Please note that some source functions (#Source-codes) are adapted from:
Bethlehem, R.A.I., Seidlitz, J., White, S.R. et al. Brain charts for the human lifespan. *Nature* 604, 525–533 (2022). [https://doi.org/10.1038/s41586-022-04554-y](https://doi.org/10.1038/s41586-022-04554-y).

This repository does not contain the original datasets. We do not have permission to distribute them. Researchers interested in accessing the full datasets should contact the authors directly.

---

## 2. Datasets

Example datasets for normative curve fitting, application, and calibration are provided for demonstration purposes only. They are intended for script testing and validation, not for research or clinical applications. For multicenter datasets used in published studies, please contact the corresponding authors for data sharing permissions.

---

## 3. Models

While we cannot share the individual-level data, we provide the outcome models (normative references and classification models).

* **#GAMLSS** contains normative models fitted on hippocampal subfield volumes.
* **#DPS** contains disease-specific classification models for disease propensity score (DPS) estimation.

---

## 4. Scripts

We have provided the main scripts for normative reference development and clinical applications.

### 4.1 Required R packages and installation

**Check-and-install-packages.R**: Automatically checks and installs the required R packages for model fitting, application, and plotting.

### 4.2 Normative curve estimate and milestone determination

**Normative-model-fit.R**: Fits normative models, estimates peak ages, and generates normative growth curves.

### 4.3 Bootstrap analysis of normative curves

**Bootstrap-normative-mdoel-fit.R**: Performs bootstrap analyses of normative curves and peak ages (default 1000 iterations).

### 4.4 Apply normative model to disease datasets

**Disease-application-normative-model.R**: Applies normative models to disease datasets to calculate deviation scores.

### 4.5 Apply normative model to an individual case

**Individual-application-normative-model.R**: Applies normative models to an individual dataset for personalized deviation scoring.

### 4.6 Calibrate normative model by new dataset

**Calibration-normative-model-using-new-dataset.R**: Calibrates normative models using new control datasets to account for site/scanner variability.

### 4.7 Statistical analyses of deviation scores across diseases

**Stastical-analysis-deviations-across-diseases.R**: Performs group-level statistical comparisons of deviation scores across diseases.

4.8 Other clinical task related scripts

ROC&DPS.R: Performs receiver operating characteristic (ROC) analysis and estimates disease propensity scores (DPS) based on deviation profiles.

cognitive-score-prediciton.R: Predicts cognitive test scores from normative deviation or centile scores using regression models.

DBS.R: Builds a support vector machine (SVM) model to classify DBS treatment outcomes (effective vs. ineffective) and evaluates classification performance using ROC curves and AUC.

HC&MCI-risk-stratificaiton.R: Implements longitudinal risk stratification models to predict whether healthy controls (HC) will convert to mild cognitive impairment (MCI), and whether MCI cases will further progress to Alzheimer’s disease (AD).

MS-Prognosis-risk-stratificaiton.R: Conducts prognosis risk stratification for multiple sclerosis patients, focusing on predicting long-term disease progression.

associations_centile_scores_cognition.R: Examines the associations between hippocampal centile scores and cognitive performance across multiple domains.---

## 5. License

MIT License

Copyright (c) [2025] [Version V1.0]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



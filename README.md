# CoxKL: Kullback-Leibler-Based Cox Models

This package is for the paper "Incorporating External Risk Information with the Cox Model under Population Heterogeneity: Applications to Trans-Ancestry Polygenic Hazard Scores".

## Purpose

Incorporating external information can improve the performance of risk discrimination in an internal small-sized cohort. However, given that external risk scores and internal individual-level data come from different populations, ignoring population heterogeneity can introduce substantial bias. Thus, we develop a Kullback-Leibler-based Cox model (CoxKL) to integrate internal individual-level time-to-event data with external risk scores derived from published prediction models. 

## Citation

Di Wang, Wen Ye, Ji Zhu, Gongjun Xu, Weijing Tang, Matthew Zawistowski, Lars G. Fritsche, and Kevin He (2023). Incorporating External Risk Information with the Cox Model under Population Heterogeneity: Applications to Trans-Ancestry Polygenic Hazard Scores. arXiv:2302.11123

## Installation

```
#Install the package, need to install the devtools packages:
install.packages("devtools")
devtools::install_github("UM-KevinHe/CoxKL")
```

## Using CoxKL

```
#Use CoxKL estimate
KL_Cox_Estimate(z, delta, time, RS_internal, eta, tol=1.0e-7)
```
- z: covariate matrix.
- delta: vector of event indicators (1 = event, 0 = censor).
- time: vector of observed times.
- RS_internal: vector of predicted risk scores for subjects in the internal data.
- eta: integration weight. Note that if eta = 0, the CoxKL model reduces to a classical Cox model using internal data only.

```
#Use cross-validation to select the optimal integration weight
kl_cox(RiskScore, eta_min, eta_max, length.out, eta_minby, df_internal, criteria, fold)
```
- RiskScore: vector of predicted risk scores for subjects in the internal data.
- eta_min, eta_max, length.out, eta_minby: for the integration weight sequence consturction, where eta_min is the start value, eta_max is the end value, length.out is length of eta sequence for each round of parameter search, and eta_minby is the minimal increment of the eta sequence. 
- df_internal: internal data (Z, status, time).
- criteria: model performance criteria for cross-validation, can be chose from "V&VH", "C-Index", "LP:log-partial likelihood", "LP: C-Index".
- fold: desired number of folds for cross-validation. 
  
## Simulation example

Examples can be performed with the following tutorial:
- download and run R file "Simulation.R"



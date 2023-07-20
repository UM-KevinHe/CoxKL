# DiscreteKL: Kullback-Leibler-Based Cox Models

This package is for the paper "Incorporating External Risk Information with the Cox Model under Population Heterogeneity: Applications to Trans-Ancestry Polygenic Hazard Scores".

## Purpose

Incorporating external information can improve the performance of risk discrimination in an internal small-sized cohort. However, given that external risk scores and internal individual-level data come from different populations, ignoring population heterogeneity can introduce substantial bias. Thus, we develop a Kullback-Leibler-based Cox model (CoxKL) to integrate internal individual-level time-to-event data with external risk scores derived from published prediction models. 

## Citation

Di Wang, Wen Ye, Ji Zhu, Gongjun Xu, Weijing Tang, Matthew Zawistowski, Lars G. Fritsche, and Kevin He (2023). Incorporating External Risk Information with the Cox Model under Population Heterogeneity: Applications to Trans-Ancestry Polygenic Hazard Scores. arXiv:2302.11123

## Installation

```
#Install the package, need to install the devtools packages:
install.packages("devtools")
devtools::install_github(diwangi/CoxKL")
```

## Simulation example

Examples can be performed with the following tutorial:
- Low dimensional setting: download and run R file "Simulation.R"



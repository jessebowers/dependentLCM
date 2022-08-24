# dependentLCM

Dependent LCM is an R package running Rcpp (a C++ interface).

Latent Class Models (LCMs) are used to cluster multivariate categorical data (e.g. group participants based on survey responses). Traditional LCMs assume a property called conditional independence. This assumption can be restrictive, leading to model misspecification and overparameterization. To combat this problem, we developed a novel Bayesian model called a Dependent Latent Class Model (DLCM), which permits conditional dependence. Compared to traditional LCMs, DLCMs are effective in applications with time series, overlapping items, and structural zeroes.

Bowers, J., & Culpepper, S. (2022). Dependent Latent Class Models (Version 1). arXiv. https://doi.org/10.48550/ARXIV.2205.08677

## Installation

Can be installed directly from R using the devtools package.

```R
library(Rcpp)
library(devtools)
devtools::install_github(repo="jessebowers/dependentLCM", ref="main")
```

## Example Usage

```R
library(dependentLCM)

# Get Data
library(pks)
data(probability, package="pks")
xdf <- probability[,c("b101", "b102", "b103", "b104", "b105", "b106", "b107", "b108", "b109", "b110", "b111", "b112", "b201", "b202", "b203", "b204", "b205", "b206", "b207", "b208", "b209", "b210", "b211", "b212")]

# Run Model
set.seed(4)
dlcm <- dependentLCM_fit(
  nitr = 6000
  , df=xdf
  , nclass=3
)
dlcm$summary <- dlcm.summary(dlcm)
```
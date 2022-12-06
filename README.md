# BayesTransitionModels
[![R-CMD-check](https://github.com/AlkemaLab/BayesTransitionModels/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/AlkemaLab/BayesTransitionModels/actions/workflows/R-CMD-check.yaml)

R package for fitting Bayesian transition models for demographic and health
indicators.

# Installation

Dependencies
- [`cmdstanr`](https://mc-stan.org/cmdstanr/). Instructions for installing `cmdstanr` are available in their [Getting started](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) guide.
- [`fpemlocal`](https://github.com/AlkemaLab/fpemlocal). For now, we need to first install [`JAGS`](https://mcmc-jags.sourceforge.io/) by downloading the latest release from their website. Note: I'm working on getting around this requirement, either by changing `fpemlocal` to not require JAGS, or move the necessary parts of `fpemlocal` over to this package. Next, install `fpemlocal` from Github:
```
   remotes::install_github("AlkemaLab/fpemlocal")
```

Install `BayesTransitionModels` from Github:
```
remotes::install_github("AlkemaLab/BayesTransitionModels")
```

This work was supported, in whole or in part, by the Bill & Melinda Gates Foundation (INV-00844). 

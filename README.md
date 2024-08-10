
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaboDynamics : developmental version

<!-- badges: start -->
<!-- badges: end -->

MetaboDynamics provides a framework of Bayesian models for robust and
easy analysis of longitudinal Metabolomics Data. It takes concentration
tables as input and provides robust estimation of mean concentrations,
functional enrichment analysis employing the KEGG database and
comparison of dynamic clusters of different biological conditions.

## Installation

You can install the development version of MetaboDynamics from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("KatjaDanielzik/MetaboDynamics")
```

## Overview

MetaboDynamics facilitates the analysis of longitudinal metabolomics
data. Common tools mostly only allow the comparison between two
timepoints or experimental conditions and are using frequentist
statistical methods.

As metabolomics data is often noisy, robust methods for the estimation
of dynamics are needed. MetaboDynamics employs a hierachical Bayesian
model. MetaboDynamics allows longitudinal analysis over multiple
timepoint employing Bayesian models.

## Workflow

For a worked example see Vignette.

<figure>
<img src="/man/figures/README-MetaboDynamics-pitch.png"
alt="Workflow" />
<figcaption aria-hidden="true">Workflow</figcaption>
</figure>

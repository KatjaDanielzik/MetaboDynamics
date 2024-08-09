
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaboDynamics

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

This package was built to faciliate the analysis of longitudinal
metabolomics data. Common tools mostly only allow the comparison between
two timepoints or experimental conditions and are using frequentist
statistical methods. As metabolomics data is often noise robust methods
for the estimation of dynamics are needed. We employ a hierachical
Bayesian model for that.

## Workflow

<figure>
<img src="/man/figures/README-metabolomics_pitch_draft.png"
alt="Workflow" />
<figcaption aria-hidden="true">Workflow</figcaption>
</figure>

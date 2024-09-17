
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaboDynamics : developmental version

<!-- badges: start -->
<!-- badges: end -->

MetaboDynamics provides a framework of Bayesian models for robust and
easy analysis of longitudinal metabolomics Data. It takes concentration
tables and KEGG IDs of metabolites as input and provides robust
estimation of mean concentrations, functional enrichment analysis
employing the KEGG database and comparison of clusters of metabolite
dynamics patterns (“dynamics clusters”) under different biological
conditions.

## Installation

You can install the development version of MetaboDynamics from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("KatjaDanielzik/MetaboDynamics")
```

or from Bioconductor with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("MetaboDynamics")
```

## Overview

MetaboDynamics facilitates the analysis of longitudinal metabolomics
data. Common tools mostly only allow the comparison between two time
points or experimental conditions and are using frequentist statistical
methods.

As metabolomics data is often noisy, robust methods for the estimation
of dynamics are needed. MetaboDynamics allows longitudinal analysis over
multiple time points and experimental conditions employing 3
probabilistic models:

1)  A hierarchical Bayesian model for the robust estimation of means at
    every time point despite varying spread between time points. Its
    outputs are A) differences between time points for every metabolite,
    and B) metabolite specific dynamics profiles that can be used for
    clustering.

2)  Over-representation analysis of KEGG-functional modules in dynamics
    clusters with a quantitative model that employs a hypergeometric
    distribution and reports probabilities of a functional module being
    over-represented in a cluster.

3)  Estimation of the distances between dynamics clusters under
    different experimental conditions with a Bayesian model that infers
    the mean pairwise Euclidean distance between two clusters. In
    combination with the comparison of metabolites that compose two
    clusters this allows to spot differences and similarities between
    experimental conditions.

## Workflow

For a worked example see Vignette (Folder /vignettes, file
MetaboDynamics.html).

![](man/figures/README-MetaboDynamics_pitch.png)

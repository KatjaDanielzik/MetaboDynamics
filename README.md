
# MetaboDynamics:
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/badge/doi-10.18129/B9.bioc.MetaboDynamics%20-yellow.svg)](https://doi.org/10.18129/B9.bioc.MetaboDynamics )
[![License: GPL](https://img.shields.io/badge/license-GPL-blue.svg)](https://cran.r-project.org/web/licenses/GPL)


MetaboDynamics provides a framework of Bayesian models for robust and
easy analysis of longitudinal metabolomics Data. It takes concentration
tables with at least three replicates and KEGG IDs of metabolites as input
and provides robust estimation of mean concentrations, functional enrichment
analysis employing the KEGG database and comparison of clusters of metabolite
dynamics patterns (“dynamics clusters”) under different biological
conditions.

## Installation

MetaboDynamics is an [R](https://cran.r-project.org/) package available
from [Bioconductor](https://www.bioconductor.org). 

URL: https://www.bioconductor.org/packages/release/bioc/html/MetaboDynamics.html

To install MetaboDynamics, start R (4.5) and enter:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetaboDynamics")
```
You can also install the development version (current bug fixes and added features
can be found in the [NEWS](https://github.com/KatjaDanielzik/MetaboDynamics/blob/main/inst/NEWS.md) 
file) of MetaboDynamics from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("KatjaDanielzik/MetaboDynamics",build_vignettes=TRUE)

```

## Overview

MetaboDynamics facilitates the analysis of longitudinal metabolomics
data e.g. from untargeted LC-MS. Common tools mostly only allow the comparison
between two time points instead of analyzing the full observed dynamics profile
of metabolite concentrations over multiple time points. Furthermore common tools
mostly only allow to compare two experimental conditions and are using 
frequentist statistical methods. As metabolomics data is often noisy, robust 
methods for the estimation of mean metabolite concentrations per time point
are needed. 

MetaboDynamics allows longitudinal analysis over
multiple time points and experimental conditions employing three
probabilistic models:

1)  A hierarchical Bayesian model for the robust estimation of means at
    every time point despite varying spread between time points.Hierarchical
    Bayesian models are known to balance between over- and underfitting, allowing
    to gain as much information from the noisy data as possible while not being overly
    confident about the estimates. Its outputs are A) differences between time 
    points for every metabolite (differential concentrations),
    and B) metabolite specific dynamics profiles that can be used for
    clustering.

2)  Over-representation analysis of KEGG-functional modules such as Amino acid
    metabolism or KEGG pathways in dynamics clusters with a quantitative model 
    that employs a hypergeometric distribution and reports probabilities of a 
    functional module or pathway being over-represented in a cluster. Can also 
    estimate under-representation of functional modules.

3)  Estimation of the dynamics similarity between metabolite clusters of 
    different experimental conditions with a Bayesian model. This model infers 
    the mean pairwise Euclidean distance of composing metabolite dynamics between 
    two clusters (i.e. every metabolite dynamics from cluster A is compared with 
    every metabolite dynamics of cluster B). In combination with the comparison of 
    metabolites that compose two clusters this allows to spot differences and similarities between
    experimental conditions. For examples clusters of metabolites with similar
    metabolite composition but different dynamics between experimental conditions.

## Workflow

For a worked example on simulated data see [Vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/MetaboDynamics/inst/doc/MetaboDynamics.html) or if package is installed:

``` r
browseVignettes("MetaboDynamics")
```

![](man/figures/README-MetaboDynamics_pitch.png)

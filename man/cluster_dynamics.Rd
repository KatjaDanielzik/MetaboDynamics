% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_dynamics.R
\name{cluster_dynamics}
\alias{cluster_dynamics}
\title{cluster dynamics profiles of metabolites}
\usage{
cluster_dynamics(
  data,
  distance = "euclidean",
  agglomeration = "ward.D2",
  minClusterSize = 1,
  deepSplit = 2
)
}
\arguments{
\item{data}{result of \code{\link[=estimates_dynamics]{estimates_dynamics()}} (list of dataframes or SummarizedExperiment object)
or a list of dataframes (one dataframe per condition, list elements must be named by condition) with columns which are named "metabolite",
"mu_mean" (mean metabolite abundance log-transformed and standarized to a mean of zero and standard deviation of one per experimental condition and metabolite) and "time.ID" (numerical, specifying the experimental time point)}

\item{distance}{distance method to be used as input for hierarchical clustering \link[stats]{dist}
can be "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"}

\item{agglomeration}{agglomerative method to be used for hierarchical clustering \link[stats]{hclust}
can be "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"}

\item{minClusterSize}{minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}}

\item{deepSplit}{rough control over sensitivity of cluster analysis. Possible values are 0:4,
the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}}
}
\value{
a list with dataframes named by experimental condition or
if data is a \link[SummarizedExperiment]{SummarizedExperiment} object clustering
results are stored in metadata under "cluster"
}
\description{
convenient wrapper function for clustering of metabolite dynamics
employing the "hybrid" method of the \link[dynamicTreeCut]{dynamicTreeCut}
package for clustering and \link[stats]{hclust} for computing of distance
matrix and hierarchical clustering needed as input for dynamicTreeCut
}
\examples{
data("longitudinalMetabolomics")
data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
  longitudinalMetabolomics$metabolite \%in\% c("ATP", "L-Alanine", "GDP")]
data <- fit_dynamics_model(
  data = data,
  scaled_measurement = "m_scaled", assay = "scaled_log",
  max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
)
data <- estimates_dynamics(
  data = data, iter = 2000,
  chains = 1, condition = "condition"
)
data <- cluster_dynamics(data)
S4Vectors::metadata(data)[["cluster"]][["A"]][["data"]]

}
\seealso{
\code{\link[=fit_dynamics_model]{fit_dynamics_model()}}, \code{\link[=estimates_dynamics]{estimates_dynamics()}}, \code{\link[=plot_cluster]{plot_cluster()}}
}

\name{cluster}
\alias{cluster}
\docType{data}
\title{
Clustering solution for data_sim
}
\description{
The estimates from the dynamic model (fit_dynamics_model) were used for
hierachical clustering with 8 clusters (8 dynamic groups were simulated).
}
\usage{data("cluster")}
\format{
  A data frame with 294 observations on the following 7 variables.
  \describe{
    \item{\code{metabolite}}{Metabolite name as in data_sim}
    \item{\code{condition}}{Experimental condition as in data_sim}
    \item{\code{mu1_mean}}{estimated mean at time point 1}
    \item{\code{mu2_mean}}{estimated mean at time point 2}
    \item{\code{mu3_mean}}{estimated mean at time point 3}
    \item{\code{mu4_mean}}{estimated mean at time point 4}
    \item{\code{cluster}}{cluster ID}
  }
}

\source{
\bold{script used to obtain clustering solution}
# get distances between vectors

dd_A <- dist(
  estimates_dynamics[["A"]][, c(
    "mu1_mean", "mu2_mean",
    "mu3_mean", "mu4_mean"
  )],
  method = "euclidean"
)

# hierachical clustering

clust <- hclust(dd_A, method = "ward.D2")
clust_cut <- cutree(clust, k = 8)

# assing cluster ID to estimates

clust_A <- estimates_dynamics[["A"]][, c(
  "metabolite", "condition", "mu1_mean", "mu2_mean",
  "mu3_mean", "mu4_mean"
)]
clust_A$cluster <- clust_cut

rm(dd_A, clust, clust_cut)

# get distances between vectors

dd_B <- dist(
  estimates_dynamics[["B"]][, c(
    "mu1_mean", "mu2_mean",
    "mu3_mean", "mu4_mean"
  )],
  method = "euclidean"
)

# hierachical clustering

clust <- hclust(dd_B, method = "ward.D2")
clust_cut <- cutree(clust, k = 8)

# assing cluster ID to estimates

clust_B <- estimates_dynamics[["B"]][, c(
  "metabolite", "condition", "mu1_mean", "mu2_mean",
  "mu3_mean", "mu4_mean"
)]
clust_B$cluster <- clust_cut

rm(dd_B, clust, clust_cut)

# get distances between vectors

dd_C <- dist(
  estimates_dynamics[["C"]][, c(
    "mu1_mean", "mu2_mean",
    "mu3_mean", "mu4_mean"
  )],
  method = "euclidean"
)

# hierachical clustering

clust <- hclust(dd_C, method = "ward.D2")
clust_cut <- cutree(clust, k = 8)

# assing cluster ID to estimates

clust_C <- estimates_dynamics[["C"]][, c(
  "metabolite", "condition", "mu1_mean", "mu2_mean",
  "mu3_mean", "mu4_mean"
)]
clust_C$cluster <- clust_cut
rm(dd_C, clust, clust_cut)

cluster <- rbind(clust_A, clust_B, clust_C)
rm(clust_A, clust_B, clust_C)
}
\examples{
data(cluster)
head(cluster)
}
\keyword{datasets}

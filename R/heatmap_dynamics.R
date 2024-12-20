#' plot bubble heatmap from the numerical fit of compare_dynamics()
#'
#' @param estimates dataframe of estimates of the mean distance
#' between clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma") obtain by function compare_dynamics()
#' @param data a dataframe or containing a column specifying the metabolite
#' names to be compared and cluster IDs (column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions.
#' to be compared or a SummarizedExperiment storing the same information in
#' metadata(data) under "cluster"
#'
#' @return a bubble heat map where the color of the bubble represents the similarity
#' of two clusters in regards to their dynamics in the color and the size the
#' uncertainty of the similarity. Big bright bubbles mean high similarity with
#' low uncertainty.
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @seealso [compare_dynamics()]
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' longitudinalMetabolomics <- compare_dynamics(
#'   data = longitudinalMetabolomics,
#'   dynamics = c("mu1_mean", "mu2_mean", "mu3_mean", "mu4_mean"),
#'   cores = 1
#' )
#' heatmap_dynamics(data = longitudinalMetabolomics)
heatmap_dynamics <- function(estimates = metadata(data)[["comparison_dynamics"]][["estimates"]],
                             data) {
  # bind objects to function
  plot <- NULL
  cluster_b <- NULL
  cluster_a <- NULL
  mu_mean <- NULL
  sigma_mean <- NULL

  # input checks
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["cluster"]]
    estimates <- metadata(data)[["comparison_dynamics"]][["estimates"]]
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!is.data.frame(estimates)) {
    stop("'estimates' must be a dataframe obtained by compare_dynamics()")
  }

  # visualization
  plot <- ggplot(estimates, aes(x = cluster_b, y = cluster_a)) +
    geom_point(aes(col = 1 / mu_mean, size = ((1 / sigma_mean)))) +
    theme_bw() +
    scale_color_viridis_c(option = "viridis") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(col = "1/estimated mean", size = "1/estimated sigma") +
    ggtitle(
      "similarity of dynamics in clusters",
      "estimated mean pairwise euclidean distance"
    )
  return(plot)
}

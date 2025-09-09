#' plot bubble heatmap from the numerical fit of compare_dynamics()
#'
#' @param estimates dataframe of estimates of the mean distance
#' between clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma") obtain by function compare_dynamics()
#' @param data a dataframe or containing a column specifying the metabolite
#' names to be compared and cluster IDs (column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions.
#' to be compared or a \link[SummarizedExperiment]{SummarizedExperiment} storing the same information in
#' metadata(data) under "cluster"
#'
#' @return a bubble heat map where the color of the bubble represents the similarity
#' of two clusters in regards to their dynamics in the color and the size the
#' uncertainty of the similarity. Big bright bubbles mean high similarity with
#' low uncertainty.
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @seealso Do calculations for comparison of dynamics between clusters [compare_dynamics()]
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition %in% c("A", "B") &
#'   longitudinalMetabolomics$metabolite %in% c("ATP", "L-Alanine", "GDP")]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- estimates_dynamics(
#'   data = data
#' )
#' data <- cluster_dynamics(data, B = 1)
#' data <- compare_dynamics(
#'   data = data,
#'   cores = 1
#' )
#' S4Vectors::metadata(data)[["comparison_dynamics"]]
#' heatmap_dynamics(data = data)
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
    # bind listelements of clustering together that contain the dataframes
    data_df <- do.call(rbind, lapply(data_df, function(l) l[["data"]]))
    estimates <- metadata(data)[["comparison_dynamics"]][["estimates"]]
  }

  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (is(data_df, "tbl")) {
    data_df <- as.data.frame(data_df)
  }
  if (!is.data.frame(estimates)) {
    stop("'estimates' must be a dataframe obtained by compare_dynamics()")
  }

  # visualization
  x <- unique(data_df[order(data_df$condition, data_df$cluster), c("condition", "cluster")])

  plot <- ggplot(estimates, aes(
    x = factor(cluster_b, levels = paste0(x$condition, "_", x$cluster)),
    y = factor(cluster_a, levels = paste0(x$condition, "_", x$cluster))
  )) +
    geom_point(aes(col = 1 / mu_mean, size = ((1 / sigma_mean)))) +
    theme_bw() +
    scale_color_viridis_c(option = "viridis") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(col = "1/estimated mean", size = "1/estimated sigma") +
    xlab("cluster_b") +
    ylab("cluster_a") +
    ggtitle(
      "similarity of dynamics in clusters",
      "estimated mean pairwise euclidean distance,
       label = condition + cluster ID"
    )
  return(plot)
}

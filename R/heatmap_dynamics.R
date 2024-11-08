#' plot bubble heatmap from the numerical fit of compare_dynamics()
#'
#' @param estimates dataframe of estimates of the mean distance
#' between clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma") obtain by function compare_dynamics()
#' @param clusters a dataframe containing the dynamics and
#' cluster IDs(column named "cluster") of clusters of similar dynamics,
#' as well as a column "condition" specifying the experimental conditions
#' to be compared.
#'
#' @return a bubble heatmap where the color of the bubble represents the similarity
#' of two clusters in regards to their dynamics in the color and the size the
#' uncertainty of the similarity. Big bright bubbles mean high similarity with
#' low uncertainty.
#' @import ggplot2
#' @seealso [compare_dynamics()]
#' @export
#'
#' @examples
#' data("cluster")
#' # fit model
#' comparison <- compare_dynamics(
#'   clusters = cluster,
#'   dynamics = c("mu1_mean", "mu2_mean", "mu3_mean", "mu4_mean"),
#'   cores = 1
#' )
#' plot <- heatmap_dynamics(estimates = comparison[["estimates"]], clusters = cluster)
#' plot

heatmap_dynamics <- function(estimates, clusters) {
  posterior <- estimates

  # bind objects to function
  plot <- NULL
  cluster_b <- NULL
  cluster_a <- NULL
  mu_mean <- NULL
  "97.5%" <- NULL
  "2.5%" <- NULL
  
  # input checks
  if (!is.data.frame(estimates)) 
    stop("'estimates' must be a dataframe obtained by compare_dynamics()")
  if (!is.data.frame(clusters)) 
    stop("'clusters' must be a dataframe")
  
  # create matrix
  # how many do we have to compare ?
  x <- unique(clusters[, c("condition", "cluster")])

  # visualization
  plot <- ggplot(posterior[posterior$parameter == "mu", ], aes(x = cluster_b, y = cluster_a)) +
    geom_point(aes(col = 1 / mu_mean, size = ((1 / (`97.5%` - `2.5%`))))) +
    theme_bw() +
    scale_color_viridis_c(option = "viridis") +
    scale_x_discrete(limits = paste0(x$condition, "_", x$cluster)) +
    scale_y_discrete(limits = paste0(x$condition, "_", x$cluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(col = "1/estimated mean", size = "1/|CrI mean|") +
    ggtitle(
      "similarity of dynamics in clusters",
      "estimated mean pairwise euclidean distance"
    )
  return(plot)
}

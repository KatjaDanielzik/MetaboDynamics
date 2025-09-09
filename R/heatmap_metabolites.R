#' plot heatmap from comparison of metabolite composition compare_metabolites()
#'
#' @param distances dataframe of Jaccard indices between clusters obtained
#' by function compare_metabolites(). If compare_metabolites() was executed on
#' as SummarizedExperiment or a \link[SummarizedExperiment]{SummarizedExperiment} than this is stored in metadata(data) under "comparison_metabolites"
#' @param data a dataframe containing the columns "metabolite" specifying the
#' metabolite names to be compared and cluster IDs(column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions
#' to be compared
#' @seealso Do calculations for comparison of metabolites between clusters [compare_metabolites()]
#' @return a heatmap where the color of the tile represents the similarity
#' of two clusters in regards to their metabolite composition. The brighter
#' the color the more similar the metabolite compositions.
#' @export
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @examples
#' data("longitudinalMetabolomics")
#' longitudinalMetabolomics <- compare_metabolites(
#'   data = longitudinalMetabolomics
#' )
#' heatmap_metabolites(data = longitudinalMetabolomics)
heatmap_metabolites <- function(distances = metadata(data)[["comparison_metabolites"]],
                                data) {
  # bind variables to function
  plot <- NULL
  Jaccard <- NULL
  cluster_a <- NULL
  cluster_b <- NULL

  # Input checks
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object obtained by 'cluster_dynamics()'")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["cluster"]]
    # combine data of list elements
    data_df <- bind_rows(lapply(data_df, function(x) {
      return(x$data)
    }))
  }
  if (is(data, "list")) {
    # combine data of list elements
    data_df <- bind_rows(lapply(data, function(x) {
      return(x$data)
    }))
  }
  # convert potential tibbles into data frame
  if (is(data_df, "tbl")) {
    data_df <- as.data.frame(data_df)
  }
  if (is(data_df, "data.frame")) {
    data_df <- data_df
  }
  if (is(data_df, "tbl")) {
    data_df <- as.data.frame(data_df)
  }
  # input checks
  if (!is.data.frame(distances)) {
    stop("'distances' must be a dataframe obtained by compare_metabolites()")
  }

  x <- unique(data_df[order(data_df$condition, data_df$cluster), c("condition", "cluster")])

  plot <- ggplot(
    distances,
    aes(
      x = factor(cluster_b, levels = paste0(x$condition, "_", x$cluster)),
      y = factor(cluster_a, levels = paste0(x$condition, "_", x$cluster)), fill = Jaccard
    )
  ) +
    geom_tile(color = "white") +
    theme_bw() +
    xlab("cluster_b") +
    ylab("cluster_a") +
    scale_fill_viridis_c(option = "viridis") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(
      "similarity of metabolites in clusters",
      "metabolites=intersection/union of metabolites in cluster,
            label = condition + cluster ID"
    )

  return(plot)
}

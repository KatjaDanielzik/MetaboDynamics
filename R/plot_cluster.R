#' visualize clustering solution of cluster_dynamics()
#'
#'
#' @param data result of [cluster_dynamics()] function
#'
#' @import ggplot2
#' @import tidyr
#' @importFrom stats prcomp
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom dendextend color_branches
#' @importFrom dendextend color_labels
#' @importFrom grDevices recordPlot
#' @importFrom graphics par
#' @importFrom graphics title
#'
#' @returns a list of plots. Per experimental condition one dendrogram,
#' colored by cluster and one visualization of PCA-analysis of the clustering solution.
#' Additionally one plot visualizing the clustered dynamics over all conditions
#'
#' @export
#'
#' @seealso [cluster_dynamics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' plot_cluster(longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A"])

plot_cluster <- function(data){
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object
         obtained by function cluster_dynamics()")
  }
  if (is(data, "list")) {
    lapply(data, function(data) {
      # check for dataframe format
      if (!inherits(data, "list")) {
        stop("'data' must be a list of lists obtained by cluster_dynamics() if it is not a SummarizedExperiment object")
      }
    })
  }

  # Data transformation
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["cluster"]]
  }
  if (is(data, "list")) {
    data_df <- data
  }

  # bind variables
  metabolite <- NULL
  cluster <- NULL
  condition <- NULL
  PC1 <- NULL
  PC2 <- NULL
  mean_log_cpc_scaled <- NULL
  time.h <- NULL

  plots <- lapply(data_df, function(data) {
    dynamics <- data[["dynamics"]]
    # dendrogram
    dendro <- as.dendrogram(data[["tree"]])
    # get color identifier based on cluster membership of metabolite
    colors_to_use <- as.numeric(unique(data[["data"]][, c("metabolite", "cluster")])$cluster)
    # order by dendrogram
    colors_to_use <- colors_to_use[order.dendrogram(dendro)]
    dendro <- dendextend::color_branches(dendro, col = colors_to_use)
    dendro <- dendextend::color_labels(dendro, col = colors_to_use)
    # get labels
    labels <- data[["tree"]]
    par(cex = 0.5)
    plot(dendro)
    par(cex = 1)
    title(main = paste0(
      "Cluster Dendrogram dynamicTreeCut, method= ",
      labels$method
    ), ylab = paste0(labels$dist.method, " distance"))
    dendrogram <- recordPlot()

    # PCA plot
    data[["data"]]$cluster <- as.factor(data[["data"]]$cluster)
    PCA <- prcomp(x = data[["data"]][, dynamics])
    # store PCA result
    temp <- cbind(PCA[["x"]], data[["data"]])
    variance <- as.data.frame(summary(PCA)[["importance"]])

    PCA_plot <- ggplot(temp[order(temp$PC1, temp$PC2), ], aes(x = PC1, y = PC2, col = cluster, group = cluster, fill = cluster)) +
      stat_ellipse(
        geom = "polygon",
        alpha = 0.2
      ) + # transparency of filling
      geom_point(aes(shape = cluster)) +
      # add enough shapes
      scale_shape_manual(values = 1:nlevels(temp$cluster)) +
      theme_bw() +
      # add axis labels with proportion of variance
      xlab(paste0("PC1 (", round(variance$PC1[2], digits = 3) * 100, "%)")) +
      ylab(paste0("PC2 (", round(variance$PC2[2], digits = 3) * 100, "%)")) +
      # add viridis scale
      scale_color_viridis_d() +
      scale_fill_viridis_d() +
      ggtitle("Principal component analysis of clustering solution", "color = cluster, points = metabolite")

    return(list(dendrogram = dendrogram, PCA_plot = PCA_plot))
  })

  # plot dynamics as lineplots
  temp <- do.call(rbind, lapply(data_df, function(l) l[["data"]]))
  dynamics <- data_df[[1]][["dynamics"]]
  temp <- temp %>% pivot_longer(cols = dynamics, names_to = "time.h", values_to = "mean_log_cpc_scaled")

  lineplot <- ggplot(temp, aes(
    x = as.factor(as.numeric(time.h)),
    y = mean_log_cpc_scaled, group = metabolite,
    col = metabolite
  )) +
    geom_line() +
    theme_bw() +
    facet_grid(cols = vars(cluster), rows = vars(condition)) +
    xlab("time point") +
    ylab("deviation from mean") +
    ylim(c(-2.5, 2.5)) +
    guides(col = "none") +
    scale_colour_viridis_d() +
    ggtitle("clustered dynamics", "line = metabolite, color = metabolite,
              panels = cluster, rows = condition")
  plots[["lineplot"]] <- lineplot

  return(plots)
}

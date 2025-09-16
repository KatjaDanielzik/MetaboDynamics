#' visualize clustering solution of cluster_dynamics()
#'
#'
#' @param data result of [cluster_dynamics()] function: either a list of data frames or a SummarizedExperiment object
#'
#' @import ggplot2
#' @import patchwork
#'
#' @returns a list of plots. Per experimental condition: 1) a 'bubbletree': a phylogram with numbers on nodes
#' indicating in how many bootstraps of the posterior estimates the same clustering
#' solution was generated, 2) cluster affiliation of metabolites, 3) dynamics of metabolites per cluster,
#' 4) patchwork of 1-3, 5) order of the tips in bubbletree: needed for matching cluster plots and ORA
#'
#' @export
#'
#' @seealso [cluster_dynamics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' plots <- plot_cluster(longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A"])
#'
#' plots[["trees"]][["A"]]
plot_cluster <- function(data) {
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object obtained by function cluster_dynamics()")
  }

  # Data transformation
  if (is(data, "SummarizedExperiment")) {
    data <- metadata(data)[["cluster"]]
  }
  if (!inherits(data, "list")) {
    stop("clustering solution obtained by cluster_dynamics must be stored in metadata of SummarizedExperiment under 'cluster'")
  }
  # binding of global variables
  label <- NULL
  metabolite <- NULL
  condition <- NULL
  cluster <- NULL
  time <- NULL

  # bubbletrees of bootstrapping
  trees <- list()
  for (i in names(data)) {
    temp <- data[[i]]
    trees[[i]] <- ggtree::ggtree(temp$mean_phylo, linetype = "solid") +
      ggtree::geom_point2(mapping = aes(subset = isTip == FALSE), size = 0.5, col = "black") +
      ggtree::geom_tippoint(size = 2, fill = "white", shape = 21) +
      ggtree::geom_tiplab(color = "black", as_ylab = TRUE, align = TRUE) +
      ggtree::layout_rectangular() +
      theme_bw(base_size = 10) +
      scale_x_continuous(labels = abs) +
      ggtree::geom_nodelab(
        geom = "text", color = "#4c4c4c", size = 2.75, hjust = -0.2,
        mapping = aes(label = label, subset = isTip == FALSE)
      )+
      ggtitle("Dendrogram","number on nodes = bootstrapps")
  }

  # plot dynamics as lineplots
  trees <- lapply(trees, ggtree::revts)
  tips <- list()
  clusterplots <- list()
  lineplots <- list()
  cluster_order <- list()
  
  for (i in names(trees)) {
    tree <- trees[[i]]
    t <- tree$data
    t <- t[order(t$y, decreasing = FALSE), ]
    tips[[i]] <- t$label[t$isTip == TRUE]
    temp <- data[[i]]$data
    temp <- temp %>% tidyr::pivot_longer(cols = -c(metabolite, condition, cluster),
                                         names_to = "time", values_to = "mean")
    temp$metabolite <- factor(temp$metabolite, levels = tips[[i]])
    temp$cluster <- as.factor(temp$cluster)
    
    cluster_order[[i]] <- unique(rev(temp[order(temp$metabolite), ]$cluster)) 
    temp$cluster <- factor(temp$cluster, levels = cluster_order[[i]])
    clusterplots[[i]] <- ggplot(temp, aes(y = metabolite, x = cluster, fill = cluster)) +
      geom_tile() +
      scale_fill_viridis_d(option="turbo") +
      guides(col = "cluster") +
      ylab("") +
      xlab("")+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      ggtitle("Cluster affiliation","dynamic tree cut")
    
    temp$time <- gsub("_mean","",temp$time)
    lineplots[[i]] <- ggplot(temp, aes(x = time, y = mean,col=cluster)) +
      geom_line(aes(group = metabolite)) +
      scale_color_viridis_d(option="turbo") +
      xlab("") +
      facet_grid(cluster~.) +
      theme_bw() +
      ylim(c(-2, 2)) +
      guides(col="none") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1))+
      geom_vline(aes(xintercept = as.factor(time)), col = "grey", linetype = "dashed") +
      geom_hline(aes(yintercept = 0), col = "grey", linetype = "dashed")+
      ggtitle("Dynamics","panel = cluster ID")
  }

  patchwork <- list()
  for (i in names(data)) {
    p <- trees[[i]] | clusterplots[[i]] | lineplots[[i]]
    patchwork[[i]] <- p + plot_annotation(
      paste0("Condition ", i))
  }

  return(list(
    trees = trees, clusterplots = clusterplots, lineplots = lineplots,
    patchwork = patchwork, cluster_order = cluster_order
  ))
}

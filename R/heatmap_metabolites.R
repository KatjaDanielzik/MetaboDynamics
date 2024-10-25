#' plot heatmap from comparison of metabolite composition compare_metabolites()
#'
#' @param distances dataframe of Jaccard indices between clusters obtained 
#' by function compare_metabolites()
#' @param clusters a dataframe containing the columns "metabolite" specifying the
#' metabolite names to be compared and cluster IDs(column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions
#' to be compared
#' @seealso [compare_metabolites()]
#' @return a heatmap where the color of the tile represents the similarity
#' of two clusters in regards to their metabolite composition. The brighter
#' the color the more similar the metabolite compositions.
#' @export
#' @import ggplot2
#' @examples
#' data("cluster")
#' comparison <- compare_metabolites(
#'   clusters = cluster
#' )
#' heatmap_metabolites(distances=comparison[["Jaccard"]],clusters=cluster)


heatmap_metabolites <- function(distances,clusters){
# bind variables to function 
  plot <- NULL
  Jaccard <- NULL
  cluster_a <- NULL
  cluster_b <- NULL
  
  x <- unique(clusters[, c("condition", "cluster")])
  
  plot <-   ggplot(
  distances[distances$Var1 < distances$Var2, ],
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
  ggtitle("similarity of metabolites in clusters", "metabolites=intersection/union of metabolites in cluster, label=condition+cluster ID")
  
  return(plot)
}

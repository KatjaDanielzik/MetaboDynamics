#' Comparison of metabolite sets between dynamics clusters of different experimental conditions
#'
#' Uses the Jaccard Index to compare metabolite names between dynamics clusters of
#' different experimental conditions
#'
#' @param clusters a dataframe containing the columns "metabolite" specifying the
#' metabolite names to be compared and cluster IDs(column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions
#' to be compared
#' @param metabolite column in "clusters" that specifies either metabolite name
#' or KEGG ID or some other identifier
#'
#' @import ggplot2
#'
#' @return a list holding a 1) dataframe of Jaccard indices between clusters
#' 2) a ggplot2 object visualizing the cluster comparison
#'
#' @export
#'
#' @examples
#' data("cluster")
#' comparison <- compare_metabolites(
#'   clusters = cluster
#' )
#' comparison[["plot_metabolite_comparison"]]
compare_metabolites <- function(clusters, metabolite = "metabolite") {
  # bind variables to function
  x <- NULL
  distances <- NULL
  a <- NULL
  b <- NULL
  temp_a <- NULL
  temp_b <- NULL
  id <- NULL
  Jaccard <- NULL
  cluster_a <- NULL
  cluster_b <- NULL

  # return object
  comparison_metabolites <- list()

  # helper function for Jaccard Index
  #' @keywords internal
  similarity <- function(a, b) {
    # intersection
    intersec <- length(intersect(a, b))
    # union=unique metabolites per set + intersection
    sim <- intersec / sum((length(a) - intersec), (length(b) - intersec), intersec)
    return(sim)
  }

  # create comparison dataframe
  # how many do we have to compare ?
  x <- unique(clusters[, c("condition", "cluster")])
  id <- expand.grid(seq(nrow(x)), seq(nrow(x)))
  distances <- cbind(id,
    cluster_a = rep(NA, nrow(id)), cluster_b = rep(NA, nrow(id)),
    Jaccard = rep(NA, nrow(id))
  )

  for (i in 1:nrow(id)) {
    # fill in only half of the matrix to save computational time
    # recover condition and cluster, condition = a[1], cluster=a[2]
    a <- paste0(x[id[i, ]$Var1, 1], "_", x[id[i, ]$Var1, 2])
    distances[i, ]$cluster_a <- a
    a <- unlist(strsplit(a, "_"))
    b <- paste0(x[id[i, ]$Var2, 1], "_", x[id[i, ]$Var2, 2])
    distances[i, ]$cluster_b <- b
    b <- unlist(strsplit(b, "_"))

    # create dataframe which combines every row from a with every row from b
    temp_a <- clusters[clusters$condition == a[1] & clusters$cluster == a[2], metabolite]
    temp_b <- clusters[clusters$condition == b[1] & clusters$cluster == b[2], metabolite]

    # cat(k)
    # calculate Jaccard index
    distances[i, ]$Jaccard <- similarity(temp_a, temp_b)
  }

  comparison_metabolites[["Jaccard"]] <- distances


  # visualize

  comparison_metabolites[["plot_metabolite_comparison"]] <-
    ggplot(
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
    ggtitle("similarity of metabolites in clusters", "metabolites=intersect/union of metabolites in cluster, label=condition+cluster ID")

  return(comparison_metabolites)
}

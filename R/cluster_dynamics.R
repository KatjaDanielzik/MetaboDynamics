#' cluster dynamics profiles of metabolites
#'
#' Convenient wrapper function for clustering of metabolite dynamics
#' employing the "hybrid" method of the \link[dynamicTreeCut]{dynamicTreeCut}
#' package for clustering and \link[stats]{hclust} for computing of distance
#' matrix and hierarchical clustering needed as input for dynamicTreeCut.
#' Provides bootstrapping of clustering solution from posterior estimates of the model.
#'
#' @param data data frame or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used used to fit dynamics model
#' @param fit model fit obtained by fit_dynamics_model(). Needed if data is not a SummarizedExperiment for which the model fit is saved in metadata[["dynamic_fit"]]
#' @param distance distance method to be used as input for hierarchical clustering \link[stats]{dist}
#' can be "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#' @param agglomeration agglomerative method to be used for hierarchical clustering \link[stats]{hclust}
#' can be "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param minClusterSize minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}
#' @param deepSplit rough control over sensitivity of cluster analysis. Possible values are 0:4,
#' the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}
#' @param B number of bootstraps
#'
#' @returns a list with dataframes named by experimental condition or
#' if data is a \link[SummarizedExperiment]{SummarizedExperiment} object clustering
#' results are stored in metadata under "cluster"
#'
#' @export
#'
#' @seealso [fit_dynamics_model()], [estimates_dynamics()], [plot_cluster()]
#'
#' @import SummarizedExperiment
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom S4Vectors metadata
#' @import tidyr
#' @import dplyr
#' @import ape 
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite %in% c("ATP", "L-Alanine", "GDP")]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- estimates_dynamics(
#'   data = data
#' )
#' data <- cluster_dynamics(data)
#' S4Vectors::metadata(data)[["cluster"]][["A"]][["data"]]
#'
cluster_dynamics <- function(data, fit, distance = "euclidean",
                             agglomeration = "ward.D2",
                             minClusterSize = 1, deepSplit = 2,
                             B = 1000) {
  
  # Input checks !!!!
  
  check for stanfit
  check for data
  check for >1 time point per metabolite and condition
  
  # # check for correct number of mu_mean for every time.ID and metabolite combination
  #   time_ids <- unique(data_df$time)
  #   metabolites <- unique(data$metabolite)
  #   for (metabolite in metabolites) {
  #     if (nrow(unique(data[data$metabolite == metabolite, ])) != length(time_ids)) {
  #       stop("dataframes in 'data' must have a mu_mean value for every time.ID per metabolite'")
  #     }
  #   }
  }

  # attach variables to function
  metabolite <- NULL
  time.ID <- NULL
  mu_mean <- NULL

  two steps
  1) on mean estimates get clustering solution
  2) bootstrapping on clustering solution for bubbletree


1) on mean estimates get clustering solution
    temp <- data %>%
      select(metabolite, time.ID, mu_mean) %>%
      pivot_wider(names_from = time.ID, values_from = mu_mean)
    # save dynamics columns
    dynamics <- colnames(temp)[-1]
    # get distance matrix from data without metabolite,KEGG and condition specification
    dist <- dist(temp[, -1], method = distance)
    # convert dist to distance matrix and assign metabolite names
    dist <- as.matrix(dist)
    colnames(dist) <- unique(data$metabolite)
    rownames(dist) <- unique(data$metabolite)
    # get hierarchical clustering result
    clust <- hclust(as.dist(dist), method = agglomeration)
    clust$dist.method <- distance
    # cut clustering results with dynamic tree cut
    cutclust <- cutreeDynamic(
      dendro = clust, # dendrogram from hierarchical clustering
      distM = as.matrix(dist), # distance matrix for hybrid method
      method = "hybrid",
      minClusterSize = minClusterSize,
      deepSplit = deepSplit
    )
    # assign clustering solution to data
    temp$cluster <- cutclust
    data <- data %>% select(-c(time.ID, mu_mean))
    temp <- left_join(temp, data, by = "metabolite")
    
# bootstrapping with phylograms (function adapted from snaketron/cellmig)
    bt <- .get_boot_tp(x = , distance = distance, method= method, B = B)



  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["cluster"]] <- result
    return(data)
  } else {
    # otherwise, return the list of fits
    return(result)
  }
}

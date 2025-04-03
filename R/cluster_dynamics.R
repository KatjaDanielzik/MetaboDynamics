#' cluster dynamics profiles of metabolites
#'
#' convenient wrapper function for clustering of metabolite dynamics
#' employing the "hybrid" method of the \link[dynamicTreeCut]{dynamicTreeCut}
#' package for clustering and \link[stats]{hclust} for computing of distance
#' matrix and hierarchical clustering needed as input for dynamicTreeCut
#'
#' @param data result of [estimates_dynamics()] (list of dataframes or SummarizedExperiment object)
#' or a list of dataframes (one dataframe per condition, list elements must be named by condition) with columns which are named "metabolite",
#' "mu_mean" (mean metabolite abundance log-transformed and standarized to a mean of zero and standard deviation of one per experimental condition and metabolite) and "time.ID" (numerical, specifying the experimental time point)
#' @param distance distance method to be used as input for hierarchical clustering \link[stats]{dist}
#' can be "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#' @param agglomeration agglomerative method to be used for hierarchical clustering \link[stats]{hclust}
#' can be "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param minClusterSize minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}
#' @param deepSplit rough control over sensitivity of cluster analysis. Possible values are 0:4,
#' the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}
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
#'   data = data, iter = 2000,
#'   chains = 1, condition = "condition"
#' )
#' data <- cluster_dynamics(data)
#' S4Vectors::metadata(data)[["cluster"]][["A"]][["data"]]
#'
cluster_dynamics <- function(data, distance = "euclidean",
                             agglomeration = "ward.D2",
                             minClusterSize = 1, deepSplit = 2) {
  message("Is your data normalized and standardized?
          We recommend normalization by log-transformation.
          Scaling and centering (mean=0, sd=1) should be metabolite and condition specific.")
  
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list of dataframes or a SummarizedExperiment object")
  }

  # Data transformation
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["estimates_dynamics"]]
  }

  # convert potential tibbles into data frame
  if (is(data, "list")) {
    lapply(data, function(data) {
      if (is(data, "tbl")) {
        data <- as.data.frame(data)
      }
      if (is(data, "data.frame")) {
        data <- data
      }
      # check for dataframe format
      if (!is.data.frame(data)) {
        stop("'data' must be a list of dataframes if it is not a SummarizedExperiment object")
      }
      # check for required column namess
      if (!all(c("metabolite", "mu_mean", "time.ID") %in% colnames(data))) {
        stop("dataframes in 'data' must contain columns named 'metabolite', 'mu_mean' and 'time.ID'")
      }
    })
    data_df <- data
  }

  # check for correct number of mu_mean for every time.ID and metabolite combination
  lapply(data_df, function(data) {
    time_ids <- unique(data$time.ID)
    metabolites <- unique(data$metabolite)
    for (metabolite in metabolites) {
      if (nrow(unique(data[data$metabolite == metabolite, ])) != length(time_ids)) {
        stop("dataframes in 'data' must have a mu_mean value for every time.ID per metabolite'")
      }
    }
  })

  # attach variables to function
  metabolite <- NULL
  time.ID <- NULL
  mu_mean <- NULL


  result <- lapply(data_df, function(data) {
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
    return(list(dynamics = dynamics, data = temp, distance_matrix = dist, tree = clust))
  })

  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["cluster"]] <- result
    return(data)
  } else {
    # otherwise, return the list of fits
    return(result)
  }
}

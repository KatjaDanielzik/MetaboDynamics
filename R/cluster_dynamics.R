#' cluster dynamics profiles of metabolites
#' 
#' convenient wrapper function for clustering of metabolite dynamics 
#' employing the "hybrid" method of the \link[dynamicTreeCut]{dynamicTreeCut}
#' package for clustering and \link[stats]{hclust} for computing of distance
#' matrix and hierarchical clustering needed as input for dynamicTreeCut
#'
#' @param data result of [estimates_dynamics()] (list of dataframes or SummarizedExperiment object)
#' or a list of dataframes (one dataframe per condition, list elements must be named by condition) with columns which are named "metabolite", 
#' "mu_mean" (mean metabolite abundance) and "time.ID" (numerical, specifying the experimental time point)
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
#'   longitudinalMetabolomics$metabolite %in% c("ATP","L-Alanine","GDP")]
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
#' 
cluster_dynamics <- function(data,distance="euclidean", 
                             agglomeration="ward.D2",
                             minClusterSize=1, deepSplit=4){
  
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["estimates_dynamics"]]
  }
  
  # convert potential tibbles into data frame
  if(is(data,"list")){
    lapply(data,function(data){
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data <- data
  }
    })
    data_df <- data
  }
  
  # attach variables to function
  metabolite <- NULL
  time.ID <- NULL
  mu_mean <- NULL
  
  
result <- lapply(data_df,function(data){
    temp <- data%>%select(metabolite,time.ID,mu_mean)%>%pivot_wider(names_from = time.ID,values_from = mu_mean)
    # get distance matrix from data without metabolite specification
    dist <- dist(temp[,-1],method=distance)
    # convert dist to distance matrix and assign metabolite names
    dist <- as.matrix(dist)
    colnames(dist) <- unique(data$metabolite)
    rownames(dist) <- unique(data$metabolite)
    # get hierarchical clustering result
    clust <- hclust(as.dist(dist), method=agglomeration)
    # cut clustering results with dynamic tree cut
    cutclust <- cutreeDynamic(dendro = clust, # dendrogram from hierarchical clustering
                              distM=as.matrix(dist), # distance matrix for hybrid method
                              method = "hybrid",
                              minClusterSize = minClusterSize,
                              deepSplit = deepSplit)
    # assign clustering solution to data
    temp$cluster <- cutclust
    data <- left_join(data,temp[,c("metabolite","cluster")],by="metabolite")
    return(list(data=data,distance_matrix=dist,tree=clust))
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
#' input checks for fit_dynamics_model
#' @keywords internal
.check_fit_dynamics_input <- function(model, data, metabolite,
                                     time, condition,
                                     scaled_measurement,
                                     counts, assay, chains, cores, adapt_delta,
                                     max_treedepth, iter, warmup){
  
  if(!model%in%c("scaled_log","raw_plus_counts")){
    stop("'model' must be either 'scaled_log' or 'raw_plus_counts'")
  }
  if(model=="scaled_log"){
    # hint user if data is standardized
    message("Are your metabolite concentrations normalized and standardized?
          We recommend normalization by log-transformation.
          Scaling and centering (mean=0, sd=1) should be metabolite and condition specific.")
  }
  if(model=="raw_plus_counts"){
    # hint user if data is standardized
    message("Are your metabolite concentrations non-normalized and non-scaled?
          The chosen 'raw_plus_counts' model expects raw metabolite concentrations.")
    
    # Input checks
    if (!is.data.frame(counts)) {
      stop("'counts' must be a dataframe if you chose model 'raw_plus_counts'.
                If you cannot provide cell counts choose model 'scaled_log'")
    }
    if (!all(c(time, condition, "counts") %in% colnames(counts))) {
      stop("'counts' must contain columns named 'time','condition', and 'counts'")
    }
  }
  
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  # check if all input variables are positive integers
  if (!all(vapply(c(iter, warmup, max_treedepth), function(x) is.numeric(x) && x > 0 && x %% 1 == 0, logical(1)))) {
    stop("'iter', 'warmup', and 'max_treedepth' must be positive integers")
  }
  if (!is.numeric(adapt_delta) | !(adapt_delta > 0 & adapt_delta < 1)) {
    stop("'adapt_delta' must be numeric and in the range (0;1)")
  }
  # check if all input variables are character vectors
  if (!all(vapply(list(metabolite, time, condition, scaled_measurement), is.character, logical(1)))) {
    stop("'metabolite', 'time', 'condition', and 'scaled_measurement' must be a character vector specifying a column name of data")
  }
}

#' euclidean distance
#' compare_dynamics()
#' @param a a numeric vector
#' @param b a numeric vector of same length as a
#' @importFrom stats dist
#' @return euclidean distance between vectors
#' @keywords internal

.eu <- function(a, b) {
  temp <- rbind(a, b)
  dist <- stats::dist(temp, method = "euclidean")
  return(dist)
}

#' Jaccard Index: intersection/union
#' compare_metabolites()
#' @param a a vector
#' @param b a vector
#'
#' @return Jaccard Index of a and b
#' @keywords internal
#'
.similarity <- function(a, b) {
  # test which vector is bigger and compare smaller to bigger
  if (length(a) < length(b)) {
    # intersection
    temp <- is.element(a, b) == TRUE
  } else {
    temp <- is.element(b, a) == TRUE
  }
  intersection <- length(temp[temp == TRUE])
  # union=unique metabolites per set + intersection
  sim <- intersection / sum(
    (length(a) - intersection),
    (length(b) - intersection), intersection
  )
  return(sim)
}

#' hierarchical clustering
#'
#' @param data estimates[["mu"]]
#' @param distance distance method
#' @param agglomeration agglomeration method of hierarchical clustering
#' @param minClusterSize minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}
#' @param deepSplit rough control over sensitivity of cluster analysis. Possible values are 0:4,
#' the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}
#' @importFrom ape as.phylo
#' @importFrom dynamicTreeCut cutreeDynamic
#' @returns list of input data including clustering solution, dendrogram, phylogram
#' @keywords internal
.hierarchical_clustering <- function(data,distance,agglomeration,minClusterSize,deepSplit){
  dist <- as.matrix(dist(data[,-c(1:3)])) # - metabolite,KEGG,condition -> only times left
  # assign metabolite names
  colnames(dist) <- data$metabolite
  rownames(dist) <- data$metabolite
  # get hierarchical clustering result
  clust <- hclust(as.dist(dist), method = agglomeration)
  clust$dist.method <- distance  
  # cut clustering results with dynamic tree cut
  cutclust <- cutreeDynamic(
    dendro = clust,
    distM = as.matrix(dist),
    method = "hybrid",
    minClusterSize = minClusterSize,
    deepSplit = deepSplit
  )
  data$cluster <- cutclust
  return(list(data=data,mean_dendro=clust, mean_phylo=ape::as.phylo(x = clust)))
}


#' get bootstrapps for clustering of dynamics vectors (cluster_dynamics function)
#'
#' @keywords internal
.get_boot_ph <- function(x, distance, agglomeration, B) {
  
  e <- x
  e <- as.data.frame(e)
  # get B samples from posterior
  samples <- sample(x= seq_len(max(e$draw)),size=min(nrow(max(e$draw)),B),replace = TRUE)
  
  boot_ph <- c()
  for(i in seq_len(length(samples))) {
    
    # hclust
    temp <- e[which(e$draw==samples[i]),-c(1,2,4)] # - draw,parameter,condition
    rownames(temp) <- temp$metabolite
    hc <- hclust(dist(as.matrix(temp[,-1]), method = distance), method = agglomeration)
    boot_ph[[i]] <- as.phylo(x = hc)
  }
  return(boot_ph)
}



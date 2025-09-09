#' input checks for fit_dynamics_model
#' @param model which model to fit. Two options are available:
#' "scaled_log": taking in normalized and scaled metabolite concentrations (see scaled measurement)
#' "raw_plus_counts": tailored for in vitro untargeted LC-MS experiments, taking in "raw"
#' (i.e. not normalized and not scaled) metabolite concentrations and cell counts.
#' This model assumes independent measurement (i.e. different wells) of cell counts
#' and metabolite concentrations. Additionally it assumes that cell counts were estimated
#' e.g. by cell counters (i.e. that cells were not counted under the microscope)
#' leading to a small uncertainty of the true cell count.
#' @param data concentration table with at least three replicate measurements per
#' metabolite. Must contain columns named "metabolite" (containing names or IDs), "time" (categorical, the same for all conditions), and "condition" or colData of a \link[SummarizedExperiment]{SummarizedExperiment} object
#' Time column needs to be sorted in ascending order
#' @param scaled_measurement column of "data" that contains the concentrations per cell,
#' centered and normalized per metabolite and experimental condition (mean=0, sd=1)
#' @param counts data frame with at least one replicate per time point and condition
#' specifying the cell counts, must contain columns "time", and "condition" equivalent
#' to the specifications of "data".
#' Must contain a column named "counts" that specifies the cell counts.
#' Model assumes that the replicates of the cell counts and metabolite concentrations
#' are independent of each other (i.e. cell counts were measured in in different
#' wells than metabolite concentrations)
#' @param assay if input is a SummarizedExperiment specify the assay that should
#' be used for input, colData has to hold the columns, "condition" and "metabolite",
#' rowData the timepoint specifications, in case of the model "scaled_log"
#' assay needs to hold scaled log-transformed metabolite concentrations
#' (mean=0,sd=1 per metabolite and experimental condition), if model
#' "raw_plus_counts" is chosen must hold the non-transformed and non-scaled metabolite concentrations
#' @param cores how many cores should be used for model fitting; this
#' parallelizes the model fitting and therefore speeds it up; default=4
#' @param chains how many Markov-Chains should be used for model fitting, use at
#' least two, default=4
#' @param adapt_delta target average acceptance probability, can be adapted if
#' divergent transitions are reported, default is 0.95
#' @param max_treedepth can be adapted if model throws warnings about hitting
#' max_treedepth, warnings are mostly efficiency not validity concerns and
#' treedepth can be raised, default=10
#' @param iter how many iterations are run, increasing might help with effective
#' sample size being to low, default=2000
#' @param warmup how many iterations the model warms up for, increasing might
#' facilitate efficiency, must be at least 25% of ITER, default=iter/4
#' @returns description error messages
#' @keywords internal
.check_fit_dynamics_input <- function(model, data,
                                      scaled_measurement,
                                      counts, assay, chains, cores, adapt_delta,
                                      max_treedepth, iter, warmup) {
  if (!model %in% c("scaled_log", "raw_plus_counts")) {
    stop("'model' must be either 'scaled_log' or 'raw_plus_counts'")
  }
  if (model == "scaled_log") {
    # hint user if data is standardized
    message("Are your metabolite concentrations normalized and standardized?
          We recommend normalization by log-transformation.
          Scaling and centering (mean=0, sd=1) should be metabolite and condition specific.")
  }
  if (model == "raw_plus_counts") {
    # hint user if data is standardized
    message("Are your metabolite concentrations non-normalized and non-scaled?
          The chosen 'raw_plus_counts' model expects raw metabolite concentrations.")

    # Input checks
    if (!is.data.frame(counts)) {
      stop("'counts' must be a dataframe if you chose model 'raw_plus_counts'.
                If you cannot provide cell counts choose model 'scaled_log'")
    }
    if (!all(c("time", "condition", "counts") %in% colnames(counts))) {
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
}

#' euclidean distance
#' compare_dynamics()
#' @param a a numeric vector
#' @param b a numeric vector of same length as a
#' @importFrom stats dist
#' @returns euclidean distance between vectors
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
#' @param data_clust element "mu" of estimates
#' @param distance distance method
#' @param agglomeration agglomeration method of hierarchical clustering
#' @param minClusterSize minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}
#' @param deepSplit rough control over sensitivity of cluster analysis. Possible values are 0:4,
#' the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}
#' @importFrom ape as.phylo
#' @importFrom dynamicTreeCut cutreeDynamic
#' @returns list of input data including clustering solution, dendrogram, phylogram
#' @keywords internal
.hierarchical_clustering <- function(data_clust, distance, agglomeration, minClusterSize, deepSplit) {
  dist <- as.matrix(dist(data_clust[, -c(1:2)])) # - metabolite,condition -> only times left
  # assign metabolite names
  colnames(dist) <- data_clust$metabolite
  rownames(dist) <- data_clust$metabolite
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
  data_clust$cluster <- cutclust
  return(list(data = data_clust, mean_dendro = clust, mean_phylo = ape::as.phylo(x = clust)))
}


#' get bootstrapps for clustering of dynamics vectors (cluster_dynamics function)
#' @param x posterior of dynamics model
#' @param distance distance measure used for hierarchical clustering
#' @param agglomeration agglomeration method used for hierarchical clustering
#' @param B number of bootstrapps
#' @keywords internal
#' @returns bootsstrapps of clustering solution
.get_boot_ph <- function(x, distance, agglomeration, B) {
  e <- x
  e <- as.data.frame(e)
  # get B samples from posterior
  samples <- sample(x = seq_len(max(e$draw)), size = min(nrow(max(e$draw)), B), replace = TRUE)

  boot_ph <- c()
  for (i in seq_len(length(samples))) {
    # hclust
    temp <- e[which(e$draw == samples[i]), -c(1, 2, 4)] # - draw,parameter,condition
    rownames(temp) <- temp$metabolite
    hc <- hclust(dist(as.matrix(temp[, -1]), method = distance), method = agglomeration)
    boot_ph[[i]] <- as.phylo(x = hc)
  }
  return(boot_ph)
}

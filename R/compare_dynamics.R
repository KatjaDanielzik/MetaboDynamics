#' Comparison of metabolite dynamics clusters under different experimental conditions
#'
#' Employs a Bayesian model that assumes a normal distribution of Euclidean
#' distances between dynamics vectors (metabolite concentrations at different time points) of two clusters that come from different
#' experimental conditions to estimate the mean distance between clusters.
#'
#' @param data a dataframe or containing a column specifying the metabolite
#' names to be compared and cluster IDs (column named "cluster") of
#' clusters of similar dynamics, as well as a column "condition" specifying
#' the experimental conditions.
#' to be compared or a \link{SummarizedExperiment} storing the same information in
#' metadata(data) under "cluster"
#' @param dynamics vector specifying the column names of dataframe clusters
#' that hold the dynamics information
#' @param cores how many cores should be used for model fitting; this
#' parallelizes the model fitting and therefore speeds it up; default=4
#'
#' @importFrom rstan sampling
#' @importFrom rstan summary
#' @importFrom utils combn
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @return a list holding a 1) the model fit
#' 2) dataframe of estimates of the mean distance
#' between #' clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma"). If data input
#' was a SummarizedExperiment results are stored in metadata(data) under
#' "comparison_dynamics"
#'
#' @seealso Visualization of estimates [heatmap_dynamics()]/
#' compare metabolite composition of clusters [compare_metabolites()]
#'
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' longitudinalMetabolomics <- compare_dynamics(
#'   data = longitudinalMetabolomics,
#'   dynamics = c("mu1_mean", "mu2_mean", "mu3_mean", "mu4_mean"),
#'   cores = 1
#' )
#' S4Vectors::metadata(longitudinalMetabolomics)[["comparison_dynamics"]]
compare_dynamics <- function(data, dynamics, cores = 4) {
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["cluster"]]
  }
  # convert potential tibbles into data frame
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!is.character(dynamics)) stop("'dynamics' must be a character vector")
  
  if (!all(c("condition", "cluster") %in% colnames(data_df))) {
    stop("'data' must contain 'condition' and 'cluster' columns")
  }
  if (!all(dynamics %in% colnames(data_df))) {
    stop("All specified 'dynamics' columns must exist in `data` dataframe")
  }

  # Extract unique condition-cluster combinations
  unique_combinations <- unique(data_df[, c("condition", "cluster")])
  n_combinations <- nrow(unique_combinations)

  # Generate all pairwise combinations (upper triangle of comparison matrix)
  comparison_pairs <- combn(seq_len(n_combinations), 2)
  n_comparisons <- ncol(comparison_pairs)

  # Calculate distances for all pairs (k)
  distances <- lapply(seq_len(n_comparisons), function(k) {
    idx_a <- comparison_pairs[1, k]
    idx_b <- comparison_pairs[2, k]

    # Extract condition and cluster info for both groups
    cond_clust_a <- unique_combinations[idx_a, ]
    cond_clust_b <- unique_combinations[idx_b, ]

    # Subset the data for the two groups
    group_a <- data_df[data_df$condition == cond_clust_a$condition &
      data_df$cluster == cond_clust_a$cluster, dynamics]
    group_b <- data_df[data_df$condition == cond_clust_b$condition &
      data_df$cluster == cond_clust_b$cluster, dynamics]

    # Compute distances
    .calculate_distances(group_a, group_b, dynamics)
  })

  # Assign names to distances
  comparison_names <- apply(comparison_pairs, 2, function(idx) {
    paste(
      paste(unique_combinations$condition[idx[1]], unique_combinations$cluster[idx[1]], sep = "_"),
      "vs",
      paste(unique_combinations$condition[idx[2]], unique_combinations$cluster[idx[2]], sep = "_")
    )
  })
  names(distances) <- comparison_names

  # Prepare distance data for Stan
  N <- lengths(distances)
  max_n <- max(N)
  y_padded <- matrix(1e6, nrow = n_comparisons, ncol = max_n)
  y_padded <- do.call(rbind, lapply(seq_along(distances), function(k) {
    c(distances[[k]], rep(1e6, max_n - N[k]))
  }))

  # Fit Bayesian model
  fit <- rstan::sampling(
    object = stanmodels$m_cluster_distance,
    data = list(
      C = n_comparisons,
      N = N,
      M = max_n,
      y = y_padded
    ),
    chains = 4,
    iter = 2000,
    warmup = 500,
    algorithm = "NUTS",
    cores = cores
  )

  # Extract posterior summaries
  posterior_mu <- as.data.frame(rstan::summary(fit, pars = "mu")$summary)
  posterior_sigma <- as.data.frame(rstan::summary(fit, pars = "sigma")$summary)

  # Combine results
  posterior <- data.frame(
    comparison = comparison_names,
    mu_mean = posterior_mu$mean,
    sigma_mean = posterior_sigma$mean
  )

  # Add cluster metadata to the posterior
  cluster_info <- do.call(rbind, strsplit(comparison_names, split = " vs "))
  posterior$cluster_a <- cluster_info[, 1]
  posterior$cluster_b <- cluster_info[, 2]

  # Return results
  result <- list(
    distances = distances,
    fit = fit,
    estimates = posterior
  )
  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["comparison_dynamics"]] <- result
    return(data)
  } else {
    # otherwise, return the result
    return(result)
  }
}

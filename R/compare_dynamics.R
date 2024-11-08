#' Comparison of metabolite dynamics clusters under different experimental conditions
#'
#' Employs a Bayesian model that assumes a normal distribution of Euclidean
#' distances between dynamics vectors (metabolite concentrations at different time points) of two clusters that come from different
#' experimental conditions to estimate the mean distance between clusters.
#'
#' @param clusters a dataframe containing the dynamics and '
#' cluster IDs(column named "cluster") of clusters of similar dynamics,
#' as well as a column "condition" specifying the experimental conditions
#' to be compared.
#' @param dynamics vector specifying the column names of dataframe clusters
#' that hold the dynamics information
#' @param cores how many cores should be used for model fitting; this
#' parallelizes the model fitting and therefore speeds it up; default=4
#'
#' @importFrom rstan sampling
#' @importFrom rstan summary
#'
#' @return a list holding a 1) the model fit
#' 2) dataframe of estimates of the mean distance
#' between #' clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma")
#'
#' @seealso Visualization of estimates [heatmap_dynamics()]/
#' compare metabolite composition of clusters [compare_metabolites()]
#'
#' @export
#'
#' @examples
#' data("cluster")
#' comparison <- compare_dynamics(
#'   clusters = cluster,
#'   dynamics = c("mu1_mean", "mu2_mean", "mu3_mean", "mu4_mean"),
#'   cores = 1
#' )
#' head(comparison[["estimates"]])

compare_dynamics <- function(clusters, dynamics, cores = 4) {
  # bind objects to function
  posterior_mu <- NULL
  posterior_sigma <- NULL
  posterior <- NULL
  id <- NULL
  temp_a <- NULL
  temp_b <- NULL
  distance <- NULL
  N <- NULL
  y_padded <- NULL
  
  # Input checks
  if (!is.data.frame(clusters)) stop("'clusters' must be a dataframe")
  if (!is.character(dynamics)) stop("'dynamics' must be a character vector")
  if (!all(c("condition", "cluster") %in% colnames(clusters))) {
    stop("'clusters' must contain 'condition' and 'cluster' columns")
  }
  if (!all(dynamics %in% colnames(clusters))) {
    stop("All specified 'dynamics' columns must exist in `clusters` dataframe")
  }
  
  # return object: stores distances, fit and extracted estimates
  comparison <- vector("list",length=3)
  names(comparison) <- c("distances","fit","estimates")
  
  # create matrix
  # how many do we have to compare?
  x <- unique(clusters[, c("condition", "cluster")])
  n_comparisons <- (nrow(x) * (nrow(x) - 1))/2  # Number of pairwise comparisons
  
  # Pre-allocate distances list
  distances <- vector("list", length = n_comparisons)
  comparison_names <- vector("character", length = n_comparisons)  # Names for distances
  index <- 1  # To track insertion index
  
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(x))) {
      # fill in only half of the matrix to save computational time
      if (j > i) {
        # recover condition and cluster, condition = a[1], cluster = a[2]
        a <- paste0(x[i, 1], "_", x[i, 2])
        a <- unlist(strsplit(a, "_"))
        b <- paste0(x[j, 1], "_", x[j, 2])
        b <- unlist(strsplit(b, "_"))
        
        # create dataframe which combines every row from a with every row from b
        temp_a <- clusters[clusters$condition == a[1] & clusters$cluster == a[2], dynamics]
        temp_b <- clusters[clusters$condition == b[1] & clusters$cluster == b[2], dynamics]
        id <- expand.grid(seq_len(nrow(temp_a)), seq_len(nrow(temp_b)))
        
        distance <- data.frame(cbind(id, euclidean = NA))
        # iterate through every combination
        for (k in seq_len(nrow(id))) {
          # calculate euclidean distance between vectors
          distance[k, "euclidean"] <- .eu(temp_a[id$Var1[k], ], temp_b[id$Var2[k], ])
        }
        
        # Store the mean Euclidean distances
        comparison_name <- paste0(
          paste0(x[i, 1], "_", x[i, 2]),
          "vs", paste0(x[j, 1], "_", x[j, 2])
        )
        
        distances[[index]] <- distance$euclidean
        comparison_names[index] <- comparison_name
        index <- index + 1
      }
    }
  }
  
  # Assign distances with their respective names
  names(distances) <- comparison_names
  comparison[["distances"]] <- distances
  
  # Get number of observations
  N <- sapply(distances, length)
  
  # Vector padding with a dummy value far from expected values
  y_padded <- matrix(1e6, nrow = length(distances), ncol = max(N))
  for (i in seq_len(length(distances))) {
    y_padded[i, 1:N[i]] <- distances[[i]]
  }
  
  # Fit model to estimate mean distance and standard deviation
  fit <- rstan::sampling(
    object = stanmodels$m_cluster_distance,
    data = list(
      C = length(distances),
      N = N,
      M = max(N),
      y = y_padded
    ),
    chains = 4,
    iter = 2000,
    warmup = 500,
    algorithm = "NUTS",
    cores = cores
  )
  comparison[["fit"]] <- fit
  
  # Prepare posterior for visualization
  posterior_mu <- as.data.frame(rstan::summary(fit, pars = c("mu"))$summary)
  posterior_mu$parameter <- "mu"
  posterior_mu$comparison <- comparison_names
  posterior_mu$mu_mean <- posterior_mu$mean
  posterior_sigma <- as.data.frame(rstan::summary(fit, pars = c("sigma"))$summary)
  posterior_sigma$parameter <- "sigma"
  posterior_sigma$comparison <- comparison_names
  posterior_sigma$sigma_mean <- posterior_sigma$mean
  posterior <- cbind(posterior_mu, sigma_mean = posterior_sigma$sigma_mean)
  posterior$cluster_a <- do.call(
    rbind,
    strsplit(posterior$comparison, split = "vs")
  )[, 1]
  posterior$cluster_b <- do.call(
    rbind,
    strsplit(posterior$comparison, split = "vs")
  )[, 2]
  
  comparison[["estimates"]] <- posterior
  
  return(comparison)
}

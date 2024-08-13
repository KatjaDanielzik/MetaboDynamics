#' Comparison of dynamic between clusters of different experimental conditions
#'
#' Employs a Bayesian model that assumes a normal distribution of euclidean
#' distances between dynamic vectors of two clusters that come from different
#' experimental conditions to estimate the mean distance between clusters.
#'
#' @param clusters a dataframe containing the dynamics and '
#' cluster IDs(column named "cluster") of clusters of similar dynamics,
#' as well as a column "condition" specifying the experimental condtions
#' to be compared.
#' @param dynamics vector specifying the column names of dataframe clusters
#' that hold the dynamic information
#'
#' @importFrom stats dist
#' @importFrom rstan sampling
#' @importFrom rstan summary
#' @import ggplot2
#'
#' @return a list holding a 1) dataframe of estimates of the mean distance
#' between #' clusters of different experimental conditions ("mean") and the
#' standard deviation ("sigma")
#' 2) the model fit
#' 3) a ggplot2 object visualizing the cluster comparison
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compare_dynamics <- functions(
#'   clusters = cluster,
#'   dynamics = c("mu1_mean", "mu2_mean", "mu3_mean", "mu4_mean"))}

compare_dynamics <- function(clusters, dynamics) {
  # bind objects to function
  posterior_mu <- NULL
  posterior_sigma <- NULL
  posterior <- NULL
  distances <- NULL
  id <- NULL
  temp_a <- NULL
  temp_b <- NULL
  distance <- NULL
  N <- NULL
  y_padded <- NULL
  cluster_b <- NULL
  cluster_a <- NULL
  mu_mean <- NULL
  "97.5%" <- NULL
  "2.5%" <- NULL


  # return object
  comparison <- list()

  # helper function for distance matrix
  #' @keywords internal
  eu <- function(a, b) {
    temp <- rbind(a, b)
    dist <- stats::dist(temp, method = "euclidean")
    return(dist)
  }

  # create matrix
  # how many do we have to compare ?
  x <- unique(clusters[, c("condition", "cluster")])

  distances <- list()
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(x)) {
      # fill in only half of the matrix to save computational time
      if (j > i) {
        # cat(i,j)
        # recover condition and cluster, condition = a[1], cluster=a[2]
        a <- paste0(x[i, 1], "_", x[i, 2])
        a <- unlist(strsplit(a, "_"))
        b <- paste0(x[j, 1], "_", x[j, 2])
        b <- unlist(strsplit(b, "_"))

        # create dataframe which combines every row from a with every row from b
        temp_a <- clusters[clusters$condition == a[1] & clusters$cluster == a[2], dynamics]
        temp_b <- clusters[clusters$condition == b[1] & clusters$cluster == b[2], dynamics]
        id <- expand.grid(seq(nrow(temp_a)), seq(nrow(temp_b)))

        distance <- data.frame(cbind(id, euclidean = NA))
        # iterate through every combination
        for (k in 1:nrow(id)) {
          # cat(k)
          # calculate euclidean distance between vectors
          distance[k, ]$euclidean <- eu(temp_a[id$Var1[k], ], temp_b[id$Var2[k], ])
        }
        # assign mean euclidean distance to comparison matrix
        distances[[paste0(
          paste0(x[i, 1], "_", x[i, 2]),
          "vs", paste0(x[j, 1], "_", x[j, 2])
        )]] <- distance$euclidean
      }
    }
  }
  rm(i, j, k, id, temp_a, temp_b, a, b, distance)

  comparison[["distances"]] <- distances

  # get number of observations
  N <- vector()
  for (i in 1:length(names(distances)))
  {
    N[i] <- length(distances[[i]])
  }
  # because we have different vector lengths -> vector padding is needed
  # use as dummy value something that is far away from expected values
  y_padded <- matrix(1e6, nrow = length(names(distances)), ncol = max(N))
  for (i in 1:length(names(distances))) {
    y_padded[i, 1:N[i]] <- distances[[i]]
  }

  # fit model to estimate mean distance and standard deviation of distances
  # between two clusters
  fit <- rstan::sampling(
    object = stanmodels$m_cluster_distance,
    data = list(
      C = length(names(distances)),
      N = N,
      M = max(N),
      y = y_padded
    ),
    chains = 4,
    iter = 2000,
    warmup = 500,
    algorithm = "NUTS",
    cores = 4
  )
  rm(y_padded, N)
  comparison[["fit"]] <- fit

  # prepare posterior for visualization
  posterior_mu <- as.data.frame(rstan::summary(fit, pars = c("mu"))$summary)
  posterior_mu$parameter <- c(rep("mu", nrow(posterior_mu)))
  posterior_mu$comparison <- names(distances)
  posterior_mu$mu_mean <- posterior_mu$mean
  posterior_sigma <- as.data.frame(rstan::summary(fit, pars = c("sigma"))$summary)
  posterior_sigma$parameter <- c(rep("sigma", nrow(posterior_sigma)))
  posterior_sigma$comparison <- names(distances)
  posterior_sigma$sigma_mean <- posterior_sigma$mean
  posterior <- cbind(posterior_mu, sigma_mean = posterior_sigma$sigma_mean)
  posterior$cluster_a <- do.call(
    rbind,
    strsplit(
      x = posterior$comparison,
      split = "vs"
    )
  )[, 1]
  posterior$cluster_b <- do.call(
    rbind,
    strsplit(
      x = posterior$comparison,
      split = "vs"
    )
  )[, 2]
  rm(fit, posterior_mu, posterior_sigma)

  comparison[["estimates"]] <- posterior

  # visualization
  comparison[["plot_dynamic_comparison"]] <-
    ggplot(posterior[posterior$parameter == "mu", ], aes(x = cluster_b, y = cluster_a)) +
    geom_point(aes(col = 1 / mu_mean, size = ((1 / (`97.5%` - `2.5%`))))) +
    theme_bw() +
    scale_color_viridis_c(option = "viridis") +
    scale_x_discrete(limits = paste0(x$condition, "_", x$cluster)) +
    scale_y_discrete(limits = paste0(x$condition, "_", x$cluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(col = "1/estimated mean", size = "1/|CrI mean|") +
    ggtitle(
      "similarity of dynamics in clusters",
      "estimated mean pairwise distance"
    )
  rm(x)
  return(comparison)
}

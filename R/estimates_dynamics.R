#' Extracts parameter estimates from numeric fit of Bayesian model of dynamics
#'
#' Extracts the mean concentrations (mu) at every timepoint from the dynamics model fit, the 95% highest density interval (HDI), the estimated standard deviation of metabolite concentrations at every time point (sigma), and the pooled standard deviation of every metabolite over all timepoints (lambda).
#' Additionally samples from the posterior of mu can be drawn. This can be helpful if p.e. one wants to estimate the clustering precision. Lambda can be used for clustering algorithms such as VSClust that also take the variance into account.
#'
#' @param data dataframe or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used used to fit dynamics model, must contain a column specifying KEGG IDs, column named "condition" specifiyng the experimental condition and a column named "time" specifying the timepoints.
#' If it is a SummarizedExperiment object the dynamic fits must be stores in metadata(data)
#' under "dynamic_fits"
#' @param assay of the SummarizedExperiment object that was used to fit the dynamics
#' model
#' @param kegg column in "data" that contains the KEGG IDs or other identifier of metabolites
#' @param condition name of column in dataframe data that specifies the experimental condition
#' @param time column in "data" that contains the time point identifiers
#' @param condition name of column in dataframe data that specifies the experimental condition
#' @param fits list of model fits for which estimates should be extracted
#' @param iter how many iterations were used to fit the dynamics model
#' @param warmup how many warm-up iterations were used to fit the dynamics model
#' @param chains how many chains were used to fit the dynamics model
#' @param samples how many posterior samples should be drawn (p.e. for check of clustering precision)
#'
#' @seealso Fit the dynamic model [fit_dynamics_model()].
#' Diagnostics of the dynamic model [diagnostics_dynamics()]
#' Visualization of estimates with [plot_estimates()]
#'
#' @import dplyr
#' @importFrom stats runif
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @return a list of dataframes (one per experimental condition) that contains
#' the estimates at the timepoints and samples from the posterior
#' (number as specified in samples), delta_mu specifies the difference between
#' time point specified in column "time.ID" and subsequent time point (delta_mu
#' in row time.ID=1: mu(time point 2)- mu(time point 1)) if number of time points
#' in dataset is >1
#' If data is a summarizedExperiment object the estimates are stored in
#' metadata(data) under "estimates_dynamics"
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite == "ATP"]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- estimates_dynamics(
#'   data = data, iter = 2000,
#'   chains = 1, condition = "condition"
#' )
#  S4Vectors::metadata(data)[["estimates_dynamics"]]
#'
estimates_dynamics <- function(data, assay = "scaled_log",
                               kegg = "KEGG", condition = "condition", time = "time",
                               fits = metadata(data)[["dynamic_fits"]],
                               iter = 2000, warmup = iter / 4, chains = 4,
                               samples = 1) {
  # Input checks
  if (!is.data.frame(data) & !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or colData of a SummarizedExperiment object")
  }
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    t <- nrow(rowData(data))
    data_df <- as.data.frame(cbind(
      as.data.frame(t(assays(data)[[assay]])),
      as.data.frame(colData(data))
    ))
    data_df <- data_df %>% pivot_longer(
      cols = seq_len(t), names_to = "time",
      values_to = "scaled_measurement"
    )
    fits <- metadata(data)[["dynamic_fits"]]
  }
  
  # convert potential tibbles into data frame
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  # check if all elements of fits are stanfit objects
  if (!all(vapply(fits, function(x) inherits(x, "stanfit"), logical(1)))) {
    stop("'fits' must be a list of stanfit objects")
  }
  if (!all(vapply(list(time, condition, kegg), is.character, logical(1)))) {
    stop("'time', 'kegg' and 'condition' must be a character vector specifying a column name of data")
  }
  if (!all(c(time, kegg, condition) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'condition', 'kegg' and 'time'")
  }
  if (!(all(c(warmup, iter, chains, samples) > 0 & c(warmup, iter, chains, samples) %% 1 == 0))) {
    stop("'iter', 'warmup', 'chains', and 'samples' must be positive integers")
  }

  M <- length(unique(data_df$metabolite))
  t <- length(unique(data_df[[time]]))
  # Unique conditions
  conditions <- unique(data[[condition]])

  # Apply helper function to all conditions
  dynamics <- lapply(conditions, function(cond) {
    fit <- fits[[cond]]
    pS <- as.data.frame(rstan::summary(fit)$summary)
    mu_posterior <- as.matrix(fit)

    # Precompute random draws
    random_draws <- floor(runif(samples, min = 1, max = (iter - warmup) * chains))

    # Generate all parameter combinations
    mu_indices <- expand.grid(metabolite = seq_len(M), timepoint = seq_len(t))
    mu_indices_names <- paste0(
      "mu[", mu_indices$metabolite,
      ",", mu_indices$timepoint, "]"
    )

    sigma_indices <- expand.grid(metabolite = seq_len(M), timepoint = seq_len(t))
    sigma_indices_names <- paste0(
      "sigma[", sigma_indices$metabolite,
      ",", sigma_indices$timepoint, "]"
    )
    lambda_indices_names <- paste0("lambda[", seq_len(M), "]")
    delta_mu_indices <- expand.grid(metabolite = seq_len(M), timepoint = seq_len(t - 1))
    delta_mu_indices_names <- paste0(
      "delta_mu[", delta_mu_indices$metabolite,
      ",", delta_mu_indices$timepoint, "]"
    )

    # Extract and format results
    mu_data <- data.frame(
      metabolite.ID = rep(seq_len(M), each = t),
      time.ID = rep(seq_len(t), M),
      condition = cond,
      mu_mean = pS[mu_indices_names, "mean"],
      mu_lower = pS[mu_indices_names, "2.5%"],
      mu_higher = pS[mu_indices_names, "97.5%"]
    )

    sigma_data <- data.frame(
      sigma_mean = pS[sigma_indices_names, "mean"],
      sigma_lower = pS[sigma_indices_names, "2.5%"],
      sigma_higher = pS[sigma_indices_names, "97.5%"]
    )

    lambda_data <- data.frame(
      metabolite.ID = seq_len(M),
      lambda_mean = pS[lambda_indices_names, "mean"],
      lambda_lower = pS[lambda_indices_names, "2.5%"],
      lambda_higher = pS[lambda_indices_names, "97.5%"]
    )

    if (t > 1) {
      delta_mu_data <- data.frame(
        metabolite.ID = rep(seq_len(M),t-1),
        time.ID = rep(seq_len(t - 1), each=M), # one less delta_t than timepoints
        delta_mu_mean = pS[delta_mu_indices_names, "mean"],
        delta_mu_lower = pS[delta_mu_indices_names, "2.5%"],
        delta_mu_higher = pS[delta_mu_indices_names, "97.5%"]
      )
    }

    # Combine results into a single data frame
    result <- cbind(mu_data, sigma_data)

    # Add lambda and delta_mu values
    result <- left_join(result, lambda_data, by = "metabolite.ID")
    if (t > 1) {
      result <- left_join(result, delta_mu_data, by = c("metabolite.ID", "time.ID"))
    }
    # Add samples
    sample_data <- vapply(random_draws, function(draw) {
      vapply(mu_indices_names, function(name) mu_posterior[draw, name], numeric(1))
    }, numeric(length(mu_indices_names)))
    colnames(sample_data) <- paste0("mu_sample_", seq_len(samples))
    rownames(sample_data) <- NULL
    result <- cbind(result, sample_data)

    # Map metabolite IDs to KEGG and names
    metabolites <- unique(data_df[c("metabolite", kegg)])
    metabolites$metabolite.ID <- as.numeric(as.factor(metabolites$metabolite))

    # Join metadata
    result <- left_join(result, metabolites, by = "metabolite.ID")
    result <- result %>% relocate("metabolite", .after = "metabolite.ID")
    result <- result %>% relocate(kegg, .after = "metabolite")
    result <- unique(result)

    return(result)
  })

  names(dynamics) <- conditions

  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["estimates_dynamics"]] <- dynamics
    return(data)
  } else {
    # otherwise, return the list of fits
    return(dynamics)
  }
}

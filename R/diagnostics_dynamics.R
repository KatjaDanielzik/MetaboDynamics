#' Extracts diagnostic criteria from numeric fit of Bayesian model of dynamics
#'
#' gathers number of divergences, rhat values, number of effective
#' samples (n_eff) and provides plots for diagnostics criteria as well as
#' posterior predictive checks. Output dataframe "model_diagnostics" contains
#' information about experimental condition, number of divergent transitions
#' and rhat and neff values for all timepoints.
#' @param data dataframe or a \link[SummarizedExperiment]{SummarizedExperiment} used to fit dynamics model
#' column of "time" that contains time must be numeric, has to contain columns
#' specifying the metabolite named "metabolite", and column specifiying the time
#' point named "time", a column named "condition" must specify the experimental condition. 
#' @param assay of the SummarizedExperiment object that was used to fit the dynamics
#' model
#' @param fit model fit for which diagnostics should be extracted, is the
#' object that gets returned by fit_dynamics_model(), or if a SummarizedExperiment
#' object the results of fit_dynamics_model() are stored in metadata(data) under
#' "dynamic_fit"
#' @param warmup number of warmup iterations used for model fit
#' @param iter number of iterations used for model fit
#' @param chains number of chains used for model fit
#'
#' @seealso [estimates_dynamics()]
#' parent function [fit_dynamics_model()]
#' visualization functions: [plot_diagnostics()]/[plot_PPC()]
#'
#' @import tidyr
#' @import dplyr
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom stats setNames
#'
#' @return a list which contains diagnostics criteria of all conditions in a
#' dataframe (named "model_diagnostics") and one dataframe per condition that
#' contains necessary information for Posterior predictive check
#' (named "PPC_condition"). If data is a summarizedExperiment object the diagnostics
#' are stored in metadata(data) "diagnostics_dynamics"
#' @export
#'
#' @examples
#'data("longitudinalMetabolomics")
#'data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition %in% c("A","B") &
#'                                   longitudinalMetabolomics$metabolite =="ATP"]
#'data <- fit_dynamics_model(
#'  data = data,
#'  scaled_measurement = "m_scaled", assay = "scaled_log",
#'  max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#')
#'data <- diagnostics_dynamics(
#'  data = data, assay = "scaled_log",
#'  iter = 2000, chains = 1,
#'  fit = metadata(data)[["dynamic_fit"]]
#')
#'S4Vectors::metadata(data)[["diagnostics_dynamics"]][["model_diagnostics"]]
#'S4Vectors::metadata(data)[["diagnostics_dynamics"]][["posterior"]]
#' 

diagnostics_dynamics <- function(data, assay = "scaled_log",
                                 iter = 2000, warmup = iter / 4, chains = 4,
                                 fit = metadata(data)[["dynamic_fit"]]) {
  # binding of variables to function
  . <- NULL
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
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
    fit <- metadata(data)[["dynamic_fit"]]
  }

  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  # read out number of metabolites and time points
  M <- length(unique(data_df$metabolite))
  N <- nrow(data_df)
  t <- length(unique(data_df$time))
  C <- length(unique(data_df$condition))
  
  # check if all elements of fits are stanfit objects
  if (!inherits(fit, "stanfit")) {
    stop("'fit' must be a stanfit object")
  }
  if (!all(c("time") %in% colnames(data_df))) {
    stop("'data' must contain a column named 'time'")
  }
  if (!(all(c(warmup, iter, chains) > 0 & c(warmup, iter, chains) %% 1 == 0))) {
    stop("'iter', 'warmup', and 'chains' must be positive integers")
  }

  # Diagnostic column names
  rhat_cols <- paste0("rhat_mu", seq_len(t), "_mean")
  neff_cols <- paste0("neff_mu", seq_len(t), "_mean")
  diag_colnames <- c(
    "metabolite.ID", "condition", "divergences",
    "treedepth_error", rhat_cols, neff_cols
  )

    # Extract model diagnostics
    divergences <- rstan::get_num_divergent(fit)
    treedepth_errors <- rstan::get_num_max_treedepth(fit)
    rhat <- rstan::summary(fit)$summary[, "Rhat"]
    n_eff <- rstan::summary(fit)$summary[, "n_eff"]

    # Create index mapping metabolites to timepoints
    metabolite_indices <- rep(seq_len(M), each = t)

    # Create data frame with all diagnostics
    diag_data <- data.frame(
      metabolite.ID = rep(seq_len(M),each=C),
      condition = rep(unique(data_df$condition),M),
      divergences = divergences,
      treedepth_error = treedepth_errors
    )
    # only extract rhat and n_eff for mu
    start <- 1
    end <- M*t*C
    diag_data <- cbind(
      diag_data,
      matrix(rhat[start:end],
        nrow = C*M, byrow = TRUE,
        dimnames = list(NULL, rhat_cols)
      ),
      matrix(n_eff[start:end],
        nrow = C*M, byrow = TRUE,
        dimnames = list(NULL, neff_cols)
      )
    )

  draws <- (iter-warmup)*chains
  # if model without cell counts y_rep else maven_rep (raw metabolite concentrations rep)
  if(fit@model_name=="m_ANOVA_partial_pooling_euclidean_distance"){
  # Posterior predictive checks for all fits
    posterior <- as.data.frame(fit, pars = "y_rep") %>%
      pivot_longer(
        cols = everything(),
        names_to = "parameter",
        values_to = "posterior"
      ) %>%
      mutate(
        metabolite.ID = rep(as.numeric(as.factor(data_df$metabolite)),draws),
        time.ID = rep(as.numeric(as.factor(data_df$time)),draws),
        condition = rep(data_df$condition,draws))
  }
  
  if(fit@model_name=="m_ANOVA_partial_pooling_cell_counts_euclidean_distance"){
    # Posterior predictive checks for all fits
    posterior <- as.data.frame(fit, pars = "maven_rep") %>%
      pivot_longer(
        cols = everything(),
        names_to = "parameter",
        values_to = "posterior"
      ) %>%
      mutate(
        metabolite.ID = rep(as.numeric(as.factor(data_df$metabolite)),draws),
        time.ID = rep(as.numeric(as.factor(data_df$time)),draws),
        condition = rep(data_df$condition,draws))
  }
  

  # Combine diagnostics and posterior predictive checks
  result <- list(model_diagnostics = diag_data,posterior = posterior)

  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["diagnostics_dynamics"]] <- result
    return(data)
  } else {
    # otherwise, return the list of fits
    return(result)
  }
}

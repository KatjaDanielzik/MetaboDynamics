#' Fits dynamic model
#'
#' MetaboDynamics provides a hierachical model for robust estimation of mean
#' concentrations over time of single metabolites at single experimental conditions.
#'
#' @param data concentration table containing the columns "metabolite", "condition", "cpc_stand" for default
#' @param metabolite column of "data" that contains the metabolite names or IDs
#' @param time column of "time" that contains time as numeric, make sure your time column is ordered from lowest to highest for the model to work
#' @param condition column of "data" that contains the experimental conditions
#' @param cpc column of "data" that contains the concentrations per cell, centered and normalized per metabolite and experimental condition (mean=0, sd=1)
#' @param cores how many cores should be used for model fitting, this parallelizes the model fitting and therefore speeds it up default=4
#' @param chains how many Markov-Chains should be used for model fitting, use at least two, default=4
#' @param adapt_delta target average acceptance probability, can be adapted if divergent transitions are reported, default is 0.95
#' @param max_treedepth can be adapted if model throws warnings about hitting max_treedepth, warnings are mostly effeciency not validity concerns and treedepth can be raised, default=10
#' @param warm_up how many iterations the model warms up for, increasing might facilitate efficiences, must be at least 25% of ITER, default=iter/4
#' @param iter how many iterations are run, increasing might help with effective sampling size being to low, default=2000
#'
#' @return returns a list of modelfits. One modelfit named fit_condition per experimental condition
#' @export
#'
#' @examples
#' fit_dynamics_model(intra)
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MetaboDynamics

fit_dynamics_model <- function(data = intra, metabolite = "metabolite", time = "time", condition = "dose", cpc = "log_cpc_stand", chains = 4, cores = 4, adapt_delta = 0.95, max_treedepth = 10, iter = 2000, warmup = iter / 4) {
  # get unique experimental conditions
  conditions <- unique(data[[condition]])
  fits <- list()
  # loop over conditions
  for (i in 1:length(conditions)) {
    cat(i)
    # subset dataframe
    temp <- data[data[[condition]] == conditions[i], ]


    # fit model
    fit <- rstan::sampling(
      object = stanmodels$m_ANOVA_partial_pooling,
      data = list(
        y = temp[[cpc]],
        t = length(unique(temp[[time]])),
        M = length(unique(temp[[metabolite]])),
        N = nrow(temp),
        # Vector of metabolite IDs
        Me = as.numeric(as.factor(temp[[metabolite]])),
        # Vector indicating which row belongs to which timestep
        X = as.numeric(as.factor(as.numeric(temp[[time]])))
      ),
      chains = chains,
      iter = iter,
      # increase warmup so that algorithm probably chooses
      # smaller step sizes for sampling the posterior
      warmup = warmup,
      algorithm = "NUTS",
      cores = cores,
      # increase adapt_delta (target average proposal acceptance probability)
      # to decrease step size
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    )
    # assign condition name to fit and store in list
    fits[[conditions[i]]] <- assign(paste0("fit", "_", conditions[i]), fit)
  }
  # return list of fits
  return(fits)
  # cleanup
  rm(i, temp, conditions, fit, fits)
}

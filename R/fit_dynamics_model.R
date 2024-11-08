#' Fits dynamics model
#'
#' Employs a hierarchical model that assumes a normal distribution of
#' standardized (mean=0, sd=1) log(cpc) (cpc = concentration per cell) values for robust estimation of mean
#' concentrations over time of single metabolites at single experimental
#' conditions.
#'
#' @param data concentration table containing the columns "metabolite",
#' "condition", and "m_scaled" by default or colData of a SummarizedExperiment
#' \linkS4class{SummarizedExperiment} object
#' @param metabolite column of "data" that contains the metabolite names or IDs
#' @param time column of "time" that contains time as numeric, make sure your
#' time column is ordered from lowest to highest for the model to work
#' @param condition column of "data" that contains the experimental conditions
#' @param scaled_measurement column of "data" that contains the concentrations per cell,
#' centered and normalized per metabolite and experimental condition (mean=0, sd=1)
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
#'
#' @seealso [diagnostics_dynamics()]/[estimates_dynamics()]
#'
#' @return returns a list of model fits. One model fit named fit_condition per
#' experimental condition
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' longitudinalMetabolomics <- as.data.frame(SummarizedExperiment::colData(longitudinalMetabolomics))
#' data <- longitudinalMetabolomics[longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite == "ATP", ]
#' fits <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", condition = "condition",
#'   max_treedepth = 14, adapt_delta = 0.999, iter = 4000, cores = 1, chains = 1
#' )
#' fits
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom SummarizedExperiment colData
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MetaboDynamics

fit_dynamics_model <- function(data, metabolite = "metabolite",
                               time = "time", condition = "dose",
                               scaled_measurement = "m_scaled", chains = 4, cores = 4,
                               adapt_delta = 0.95, max_treedepth = 10,
                               iter = 2000, warmup = iter / 4) {

  # Input checks
  if (!is.data.frame(data)|inherits(data,"SummarizedExperiment")) 
    stop("'data' must be a dataframe or colData of a SummarizedExperiment object")
  if(!sapply(fits, function(x) inherits(x, "stanfit")))
    stop("'fits' must be a list of stanfit objects")
  if (!is.character(all(c(metabolite,time,condition,scaled_measurement)))) 
    stop("'metabolite', 'time', 'condition', and 'scaled_measurement' must be a character vector specifying a column name of data")
  if (!all(c(metabolite,time,condition,scaled_measurement) %in% colnames(data))) {
    stop("'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'")
  }
  if (!is.integer(c(warmup,iter,max_treedepth)|!c(warmup,iter,max_treedepth)>0)) {
    stop("'iter', 'warmup', and 'max_treedepth' must be positive integers")
  }
  if (!is.numeric(adapt_delta)|!(adapt_delta>0&adapt_delta<1)) {
    stop("adapt_delta must be numeric and in the range [0;1]")
  }
  
  
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }

  # get unique experimental conditions
  conditions <- unique(data[[condition]])
  fits <- list()
  # loop over conditions
  for (i in seq_len(length(conditions))) {
    # subset dataframe
    temp <- data[data[[condition]] == conditions[i], ]


    # fit model
    fit <- rstan::sampling(
      object = stanmodels$m_ANOVA_partial_pooling,
      data = list(
        y = temp[[scaled_measurement]],
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
}

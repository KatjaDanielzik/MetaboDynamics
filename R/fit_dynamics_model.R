#' Fits dynamics model
#'
#' Employs a hierarchical model that assumes a normal distribution of
#' standardized (mean=0, sd=1) log(cpc) (cpc = concentration per cell) values for robust estimation of mean
#' concentrations over time of single metabolites at single experimental
#' conditions.
#'
#' @param data concentration table with at least three replicate measurements per
#' metabolites containing the columns "metabolite",
#' "condition", and "m_scaled" by default or colData of a \link[SummarizedExperiment]{SummarizedExperiment} object
#' @param metabolite column of "data" that contains the metabolite names or IDs
#' @param time column of "time" that contains time as numeric, make sure your
#' time column is ordered from lowest to highest for the model to work
#' @param condition column of "data" that contains the experimental conditions
#' @param scaled_measurement column of "data" that contains the concentrations per cell,
#' centered and normalized per metabolite and experimental condition (mean=0, sd=1)
#' @param assay if input is a summarizedExperiment specify the assay that should
#' be used for input, colData has to hold the columns, "condition" and "metabolite",
#' rowData the timepoint specifications
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
#' @seealso Example data set[longitudinalMetabolomics].
#' Get model diagnostics [diagnostics_dynamics()]
#' Get model estimates [estimates_dynamics()]
#'
#' @return returns a list of model fits. One model fit named "condition" per
#' experimental condition. If input is a summarizedExperiment object the dynamic
#' fits are stored metadata(data) under "dynamic_fits"
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
#' S4Vectors::metadata(data)[["dynamic_fits"]]
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MetaboDynamics

fit_dynamics_model <- function(data, metabolite = "metabolite",
                               time = "time", condition = "condition",
                               scaled_measurement = "m_scaled",
                               assay = "scaled_log",
                               chains = 4, cores = 4,
                               adapt_delta = 0.95, max_treedepth = 10,
                               iter = 2000, warmup = iter / 4) {
  # hint user if data is standardized
  message("Is your data normalized and standardized?
          We recommend normalization by log-transformation.
          Scaling and centering (mean=0, sd=1) should be metabolite and condition specific.")
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
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    t <- nrow(rowData(data))
    data_df <- as.data.frame(cbind(
      as.data.frame(t(assays(data)[[assay]])),
      as.data.frame(colData(data))
    ))
    data_df <- data_df %>% pivot_longer(
      cols = seq_len(t), names_to = time,
      values_to = scaled_measurement
    )
  }
  
  # convert potential tibbles into data frame
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  # assign number of time points
  t <- length(unique(data_df[[time]]))

  # check if all input variables are character vectors
  if (!all(vapply(list(metabolite, time, condition, scaled_measurement), is.character, logical(1)))) {
    stop("'metabolite', 'time', 'condition', and 'scaled_measurement' must be a character vector specifying a column name of data")
  }
  if (!all(c(metabolite, time, condition, scaled_measurement) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'")
  }
  
  
  # validate at least triplicate measurements
  # count replicates per metabolite, time and condition
  grouped_data <- data_df %>%
    group_by(metabolite, time, condition) %>%
    summarise(count = n())
  if (any(grouped_data$count < 3) == TRUE) {
    stop("Input must contain at least three measurements per metabolite,
      time and experimental condition.")
  }

  # split data into a list of data frames, where each data frame corresponds to
  # a unique experimental condition
  data_split <- split(data_df, data_df[[condition]])

  # define a function to fit the model to a single data frame
  # apply the model fitting function to each data frame in the list
  fits <- lapply(data_split, function(temp) {
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
        X = as.numeric(as.factor(temp[[time]]))
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
  })

  # assign names to the model fits based on the experimental condition
  names(fits) <- names(data_split)

  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["dynamic_fits"]] <- fits
    return(data)
  } else {
    # otherwise, return the list of fits
    return(fits)
  }
}

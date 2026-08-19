#' Fits dynamics model
#'
#' Employs a hierarchical model that assumes a normal distribution of
#' standardized (mean=0, sd=1) log(cpc) (cpc = normalized metabolite abundance)
#' values for robust estimation of mean
#' concentrations over time of single metabolites at single experimental
#' conditions.
#' At least three replicates for metabolite concentrations per time point and condition are needed.
#' If cell counts are provided at least one replicate per time point and condition is needed.
#'
#' @param model which model to fit. Two options are available:
#' "scaled_log": taking in normalized and scaled metabolite concentrations (see scaled measurement)
#' "raw_plus_counts": tailored for in vitro untargeted LC-MS experiments, taking in "raw"
#' (i.e. not normalized and not scaled) metabolite concentrations and cell counts.
#' This model assumes independent measurement (i.e. different wells) of cell counts
#' and metabolite concentrations. Additionally it assumes that cell counts were estimated
#' e.g. by cell counters (i.e. that cells were not counted under the microscope)
#' leading to a small uncertainty of the true cell count.
#' @param model_option for each model ("scaled_log" and "raw_plus_counts") two
#' versions are available: "sd_per_time_point" and "sd_per_condition".
#' "sd_per_time_point" is suited for more replicates (at least triplicates) and
#' estimates standard deviations per time point and therefore allows for more
#' sensitive discrimination of changing abundances.
# "sd_per_condition" can be suitable for down to duplicates (but may brake for
#' high deviations between replicates) or for many metabolites, time points or conditions
#' as it reduces run time due to less parameters. For model details see vignette (browseVignettes("MetaboDynamics")).
#' @param data concentration table with at least three (model_option "sd_per_time_point) or two (model_option "sd_per_condition")
#' replicate measurements per
#' metabolite. Must contain columns named "metabolite" (containing names or IDs), "time" (categorical, the same for all conditions), and "condition" or colData of a \link[SummarizedExperiment]{SummarizedExperiment} object
#' Time column needs to be sorted in ascending order
#' @param scaled_measurement column of "data" that contains the abundances per cell (model "raw_plus_counts"),
#' or in case of model "scaled_log" the log-transformed metabolite abundances
#' centered and normalized per metabolite and experimental condition  to mean=0 and sd=1,
#' must be numeric
#' @param counts only model "raw_plus_counts", data frame with at least one replicate 
#' per time point and condition
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
#' @param prior_mean_abundance only for model 'raw_plus_counts', has to be >0, mu and sigma of
#' normal distribution covering the range for log-transformed metabolite abundances,
#' for details see Vignette "Prior setting".
#' @param prior_sd_abundance only for model 'raw_plus_counts', has to be >0, expected value
#' (most abundant) value of metabolite specific standard deviation of log-transformed metabolite abundances,
#' for details see Vignette "Prior setting". 
#' @param prior_counts only for model 'raw_plus_counts', has to be >0, expected value of cell counts,
#' for details see Vignette "Prior setting".
#' 
#'
#' @seealso Example data set[longitudinalMetabolomics].
#' Get model diagnostics [diagnostics_dynamics()]
#' Get model estimates [estimates_dynamics()]
#'
#' @return returns a list of model fits. One model fit named "condition" per
#' experimental condition. If input is a summarizedExperiment object the dynamic
#' fits are stored metadata(data) under "dynamic_fits". 'model' and 'model_option'
#' are also stored in metadata(data)
#' @export
#'
#' @examples
#' ## on scaled log-transformed metabolite concentrations
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition %in% c("A", "B") &
#'   longitudinalMetabolomics$metabolite == "ATP"]
#' data <- fit_dynamics_model(
#'   model = "scaled_log",
#'   data = data,
#'   assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' S4Vectors::metadata(data)[["dynamic_fit"]]
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom rlang as_name
#' @useDynLib MetaboDynamics

fit_dynamics_model <- function(model = "scaled_log",
                               model_option = "sd_per_condition",
                               data,
                               scaled_measurement = "m_scaled",
                               counts = NULL,
                               assay = "scaled_log",
                               chains = 4, cores = 4,
                               adapt_delta = 0.95, 
                               max_treedepth = 10,
                               iter = 2000, 
                               warmup = iter / 4,
                               prior_mean_abundance = c(12,5),
                               prior_sd_abundance = 2,
                               prior_counts = 1e7) {
  
  .check_fit_dynamics_input(
    model = model, model_option = model_option, data = data,
    scaled_measurement = scaled_measurement,
    counts = counts, assay = assay, chains = chains,
    cores = cores, adapt_delta = adapt_delta,
    max_treedepth = max_treedepth, iter = iter, warmup = warmup, prior_mean_abundance,
    prior_sd_abundance, prior_counts
  )

  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    t <- nrow(rowData(data))
    data_df <- as.data.frame(cbind(
      as.data.frame(t(assays(data)[[assay]])),
      as.data.frame(colData(data))
    ))
    data_df <- data_df %>% pivot_longer(
      cols = seq_len(t), names_to = "time",
      values_to = scaled_measurement
    )
  }
  # convert potential tibbles into data frame
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (is(data, "tbl")) {
    data_df <- as.data.frame(data)
  }

  if (!all(c("metabolite", "time", "condition", scaled_measurement) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'")
  }

  if (!is.numeric(data_df[[scaled_measurement]])) {
    stop("'scaled_measurement' must be numeric")
  }
  if (all(is.na(data_df[[scaled_measurement]]))) {
    stop("'scaled_measurement' cannot contain NAs")
  }

  if (model == "raw_plus_counts") {
    if (is(counts, "tbl")) {
      counts <- as.data.frame(counts)
    }
  }

  # validate at least triplicate measurements
  # count replicates per metabolite, time and condition
  if (model_option == "sd_per_time_point") {
    grouped_data <- data_df %>%
      group_by(metabolite, time, condition) %>%
      summarise(count = n())
    if (any(grouped_data$count < 3) == TRUE) {
      stop("Input must contain at least three replicates per metabolite, time point and experimental condition.")
    }
  }

  if (model_option == "sd_per_condition") {
    grouped_data <- data_df %>%
      group_by(metabolite, time, condition) %>%
      summarise(count = n())
    if (any(grouped_data$count < 2) == TRUE) {
      stop("Input must contain at least two replicates per metabolite,
      time point and experimental condition. Check diagnostics and PPC carefully
      before using estimates!")
    }
  }


  # check if all conditions and time points have cell counts
  if (model == "raw_plus_counts") {
    if (!setequal(unique(data_df$time), unique(counts$time))) {
      stop("data and counts must have the same time points")
    }
    if (!setequal(unique(data_df$condition), unique(counts$condition))) {
      stop("data and counts must have the same conditions")
    }
  }


  # Binding of global variables
  time <- NULL
  metabolite <- NULL
  condition <- NULL

  if (model == "scaled_log") {
    if (model_option == "sd_per_time_point") {
      stanmodel <- stanmodels$m_ANOVA_partial_pooling_euclidean_distance
    }
    if (model_option == "sd_per_condition") {
      stanmodel <- stanmodels$m_ANOVA_partial_pooling_euclidean_distance_robust
    }
    # fit model
    fit <- rstan::sampling(
      object = stanmodel,
      data = list(
        N = nrow(data_df),
        M = length(unique(data_df$metabolite)),
        t = length(unique(data_df$time)),
        d = length(unique(data_df$condition)),
        y = data_df[[scaled_measurement]],
        # Vector of metabolite IDs
        Me = as.numeric(as.factor(data_df$metabolite)),
        # Vector indicating which row belongs to which timestep
        X = as.numeric(as.factor(data_df$time)),
        # Condition indicator
        Do = as.numeric(as.factor(data_df$condition))
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
  }

  if (model == "raw_plus_counts") {
    if (model_option == "sd_per_time_point") {
      stanmodel <- stanmodels$m_ANOVA_partial_pooling_cell_counts_euclidean_distance
    }
    if (model_option == "sd_per_condition") {
      stanmodel <- stanmodels$m_ANOVA_partial_pooling_cell_counts_euclidean_distance_robust
    }
    # fit model
    fit <- rstan::sampling(
      object = stanmodel,
      data = list(
        N = nrow(data_df),
        M = length(unique(data_df$metabolite)),
        t = length(unique(data_df$time)),
        D = length(unique(data_df$condition)),
        maven = data_df[[scaled_measurement]],
        # Vector of metabolite IDs
        Me = as.numeric(as.factor(data_df$metabolite)),
        # Vector indicating which row belongs to which timestep
        X = as.numeric(as.factor(data_df$time)),
        # Condition indicator
        Do = as.numeric(as.factor(data_df$condition)),
        Nc = nrow(counts),
        Cc = as.numeric(counts$counts),
        X_c = as.numeric(as.factor(counts$time)),
        Do_c = as.numeric(as.factor(counts$condition),
        prior_mean_abundance = prior_mean_abundance,
        prior_sd_abundance = prior_sd_abundance,
        prior_counts = prior_counts)
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
  }


  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["dynamic_fit"]] <- fit
    metadata(data)[["model"]] <- model
    metadata(data)[["model_option"]] <- model_option
    # store priors
    if(model=="raw_plus_counts"){
      metadata(data)[["priors"]] <- list(prior_mean_abundance=prior_mean_abundance,
                                         prior_sd_abundance=prior_sd_abundance,
                                         prior_counts=prior_counts)
    }
    return(data)
  } else {
    # otherwise, return the list of fits
    return(fit)
  }
}

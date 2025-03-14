#' Plots posterior predictive check of numerical fit of Bayesian dynamics model
#'
#' @param posterior a list of one dataframe per condition that
#' contains necessary information for Posterior predictive check
#' obtained by function diagnostics_dynamics()(named "PPC_condition")
#' @param data dataframe or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used to fit dynamics model
#' @param assay of the \link[SummarizedExperiment]{SummarizedExperiment}  object that was used to fit the dynamics
#' model
#' @param scaled_measurement column name of concentration values used to model fit, should be normalized by
#' experimental condition and metabolite to mean of zero and standard deviation
#' of one
#' @return a list of visual posterior predictive check, one per experimental condition
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @export
#' @seealso  parent function [diagnostics_dynamics()]
#' visualization function for diagnostics [plot_diagnostics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite %in% c("ATP", "ADP")]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- diagnostics_dynamics(
#'   data = data, assay = "scaled_log",
#'   iter = 2000, chains = 1,
#'   fits = metadata(data)[["dynamic_fits"]]
#' )
#' plot_PPC(
#'   data = data, assay = "scaled_log"
#' )
plot_PPC <- function(
    posterior = metadata(data)[["diagnostics_dynamics"]],
    data, assay = "scaled_log",
    scaled_measurement = "scaled_measurement") {
  # bind variables to function
  time.ID <- NULL
  plots <- NULL

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
    posterior <- metadata(data)[["diagnostics_dynamics"]]
    # only select posteriors
    posterior <- posterior[-1]
  }
  # convert potential tibbles into data frame
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!inherits(posterior, "list")) {
    stop("'posterior' must be a list obtained by diagnostics_dynamics()")
  }
  if (!is.character(scaled_measurement)) {
    stop("'scaled_measurement' must be a character vector specifying a column name of data")
  }
  if (!all(c(scaled_measurement, "metabolite") %in% colnames(data_df))) {
    stop("'data' must contain a column named 'scaled_measurement' and 'metabolite'")
  }


  # prepare data for PPC
  # assign metabolite and time id to data
  PPC <- data_df
  PPC$metabolite.ID <- as.numeric(as.factor(PPC$metabolite))
  PPC$time.ID <- as.numeric(as.factor(PPC$time))
  # make scaled_measurement useable with tidy evaluaions
  scaled_measurement <- as.symbol(scaled_measurement)
  scaled_measurement <- enquo(scaled_measurement)


  # plot for every experimental condition
  plots <- list()
  for (i in names(posterior)) {
    plots[[i]] <-
      ggplot(posterior[[i]], aes(x = as.factor(time.ID))) +
      geom_violin(aes(y = posterior, x = as.factor(time.ID)), scale = "count") +
      geom_jitter(data = PPC, aes(x = as.factor(time.ID), y = !!scaled_measurement), width = 0.05) + # aes_string allows us to use predefined variables
      theme_bw() +
      ylim(-5, 5) + # we standardized data so we are not expecting much smaller or bigger values
      xlab("time point") +
      ggtitle(
        paste0(
          "Posterior predicitve check ",
          ": points within violins?"
        ),
        "violins=posterior, points=data"
      )
  }
  return(plots)
}

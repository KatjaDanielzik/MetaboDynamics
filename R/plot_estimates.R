#' Visualization of parameter estimates from numeric fit of Bayesian model of dynamics
#'
#' @param estimates a list of dataframes (one per experimental condition) that contains
#' the estimates at the timepoints and samples from the posterior generated
#' by estimates_dynamics() or if data is a \link[SummarizedExperiment]{SummarizedExperiment}  estimates
#' must be stored in metadata(data) under "estimates_dynamics"
#' @param data dataframe or \link[SummarizedExperiment]{SummarizedExperiment}  used used to fit dynamics model and extract the estimates
#' @param assay of the \link[SummarizedExperiment]{SummarizedExperiment}  object that was used to fit the dynamics
#' model
#' @param time column in "data" that contains the time point identifiers
#' @param delta_t should differences between timepoints be plotted?
#' @param dynamics should dynamics be plotted?
#' @return Visualization of differences between timepoints(delta_t) and dynamics
#' profiles of single metabolites
#' @export
#' @import ggplot2
#' @import dplyr
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @seealso parent function [estimates_dynamics()]
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
#' plot_estimates(data = data, delta_t = FALSE)
#' plot_estimates(data = data, dynamics = FALSE)
plot_estimates <- function(data,
                           estimates = metadata(data)[["estimates_dynamics"]],
                           assay = "scaled_log", time = "time",
                           delta_t = TRUE, dynamics = TRUE) {
  # bind variables to function
  condition <- NULL
  metabolite.ID <- NULL
  time.ID <- NULL
  delta_mu_mean <- NULL
  delta_mu_lower <- NULL
  delta_mu_higher <- NULL
  metabolite <- NULL
  mu_mean <- NULL
  mu_higher <- NULL
  mu_lower <- NULL
  timepoints <- NULL

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
    estimates <- metadata(data)[["estimates_dynamics"]]
  }

  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!all(vapply(list(time), is.character, logical(1)))) {
    stop("'time', 'kegg' and 'condition' must be a character vector specifying a column name of data")
  }
  if (!all(c(time) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'time'")
  }
  t <- length(unique(data_df$time))

  if (!is.list(estimates)) {
    stop("'estimates' must be a list of dataframes obtained by estimates_dynamics()")
  }
  if (!is.logical(delta_t)) {
    stop("'delta_t' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(dynamics)) {
    stop("'dynamics' must be either 'TRUE' or 'FALSE'")
  }

  plots <- list()

  # visualize
  # bind if multiple conditions are analyzed
  temp <- do.call(rbind, estimates)

  # differences between timepoints
  if (delta_t == TRUE) {
    if (t < 2) {
      stop("differences between timepoints can only be plotted if
                 dataset contains more than two time points")
    }
    temp_t <- temp %>% select(
      metabolite, condition, metabolite.ID, time.ID, delta_mu_mean, delta_mu_lower,
      delta_mu_higher
    )
    temp_t <- temp_t %>% mutate(col = ifelse(delta_mu_higher < 0, "HDI>0",
      ifelse(delta_mu_lower > 0, "HDI<0", "0inHDI")
    ))

    plots[["plot_timepoint_differences"]] <-
      ggplot(temp_t[which(!is.na(temp_t$delta_mu_mean)), ], aes(y = delta_mu_mean, x = metabolite, col = col)) +
      geom_point() +
      geom_errorbar(aes(ymin = delta_mu_lower, ymax = delta_mu_higher)) +
      ylab("delta") +
      scale_color_manual(
        values = c("black", "green", "red"),
        labels = c("0inCrI", "CrI>0", "CrI<0"), name = ""
      ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(time.ID), cols = vars(condition)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      ggtitle(
        "differences between timepoints",
        "point = mean, errorbar = 95% highest density interval (CrI),
       1 = time point 2 - time point 1"
      )
  }

  # dynamics
  if (dynamics == TRUE) {
    temp_d <- temp %>% select(metabolite, condition, metabolite.ID, time.ID, mu_mean)
    plots[["plot_dynamics"]] <-
      ggplot(temp_d, aes(
        x = as.factor(as.numeric(as.factor(time.ID))),
        y = mu_mean, group = metabolite.ID, col = metabolite
      )) +
      geom_line() +
      xlab("time point") +
      scale_color_viridis_d()+
      ylab("estimated deviation from mean concentration") +
      theme_bw() +
      theme(legend.position = "none") +
      facet_grid(rows = vars(condition)) +
      ggtitle("dynamics", "color=metabolite, row labels = condition")
  }
  return(plots)
}

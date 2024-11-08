#' Visualization of parameter estimates from numeric fit of Bayesian model of dynamics
#'
#' @param estimates a list of dataframes (one per experimental condition) that contains
#' the estimates at the timepoints and samples from the posterior generated
#' by estimates_dynamics()
#' @param data dataframe or colData of a SummarizedExperiment used used to fit dynamics model
#' @param delta_t should differences between timepoints be plotted?
#' @param dynamics should dynamics be plotted?
#' @return Visualization of differences between timepoints(delta_t) and dynamics
#' profiles of single metabolites
#' @export
#' @import ggplot2
#' @seealso [estimates_dynamics()]
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
#' estimates <- estimates_dynamics(
#'   data = data, fits = fits, iter = 4000,
#'   chains = 1, condition = "condition"
#' )
#' plot_estimates(estimates = estimates, data = data, delta_t = TRUE)
#' plot_estimates(estimates = estimates, data = data, dynamics = TRUE)
plot_estimates <- function(estimates, data, delta_t = TRUE, dynamics = TRUE) {
  
  # bind variables to function
  condition <- NULL
  metabolite.ID <- NULL
  metabolite <- NULL
  higher <- NULL
  lower <- NULL
  timepoints <- NULL
  t <- length(unique(data$time))
  plots <- list()

  # input checks
  if (!is.data.frame(data)|inherits(data,"SummarizedExperiment")) 
    stop("'data' must be a dataframe or colData of a SummarizedExperiment object")
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }
  if (!is.data.frame(estimates)) 
    stop("'estimates' must be a dataframe obtained by estimates_dynamics()")
  if (!is.logical(delta_t))
    stop("'delta_t' must be either 'TRUE' or 'FALSE'")
  if (!is.logical(dynamics))
    stop("'dynamics' must be either 'TRUE' or 'FALSE'")
  
  
  # visualize
  # bind if multiple conditions are analyzed
  temp <- estimates[[1]]
  if (length(names(estimates)) > 1) {
    for (i in 2:length(names(estimates))) {
      temp <- rbind(temp, estimates[[i]])
    }
  }

  # differences between timepoints
  if (delta_t == TRUE) {
    temp_t <- temp[, c(1:3, (6 * t + 7):(6 * t + 15))]
    temp_t <- temp_t %>% pivot_longer(
      cols = -c(condition, metabolite.ID, metabolite),
      names_to = c("timepoints", ".value"), names_sep = "_"
    )
    temp_t <- temp_t %>% mutate(col = ifelse(higher < 0, "HDI>0",
      ifelse(lower > 0, "HDI<0", "0inHDI")
    ))
    plots[["plot_timepoint_differences"]] <-
      ggplot(temp_t, aes(y = as.numeric(mean), x = metabolite, col = col)) +
      geom_point() +
      geom_errorbar(aes(ymin = lower, ymax = higher)) +
      ylab("delta") +
      scale_color_manual(
        values = c("black", "green", "red"),
        labels = c("0inCrI", "CrI>0", "CrI<0"), name = ""
      ) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_grid(rows = vars(timepoints), cols = vars(condition)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
      ggtitle("differences between timepoints")
  }

  # dynamics
  if (dynamics == TRUE) {
    temp_d <- temp[, c(1:3, 4:(t + 3))]
    temp_d <- temp_d %>% pivot_longer(
      cols = -c(condition, metabolite.ID, metabolite),
      names_to = c("timepoints", ".value"), names_sep = "_"
    )
    plots[["plot_dynamics"]] <-
      ggplot(temp_d, aes(
        x = as.factor(as.numeric(as.factor(timepoints))),
        y = mean, group = metabolite.ID, col = metabolite
      )) +
      geom_line() +
      xlab("timepoint") +
      ylab("estimated mean concentration") +
      theme_bw() +
      theme(legend.position = "none") +
      facet_grid(rows = vars(condition)) +
      ggtitle("dynamics", "color=metabolite")
  }
  return(plots)
}

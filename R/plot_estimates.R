#' Visualization of parameter estimates from numeric fit of Bayesian model of dynamics
#'
#' @param estimates a list of data frames (elements: mu, sigma, lambda, euclidean_distance) that contains
#' the model estimates by estimates_dynamics() or if data is a \link[SummarizedExperiment]{SummarizedExperiment}  estimates
#' must be stored in metadata(data) under "estimates_dynamics"
#' @param data \link[SummarizedExperiment]{SummarizedExperiment} used to fit dynamics model and extract the estimates
#' @param delta_t should differences between time points be plotted?
#' @param dynamics should dynamics be plotted?
#' @param distance_conditions should differences in metabolite specific dynamic should be plotted?
#' @return Visualization of differences between time points(delta_t) and dynamics
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
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition %in% c("A", "B") &
#'   longitudinalMetabolomics$metabolite %in% c("ATP")]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- estimates_dynamics(
#'   data = data
#' )
#' plot_estimates(data = data, delta_t = TRUE, dynamic = FALSE, distance_conditions = FALSE)
#' plot_estimates(data = data, delta_t = FALSE, dynamic = TRUE, distance_conditions = FALSE)
#' plot_estimates(data = data, delta_t = FALSE, dynamic = FALSE, distance_conditions = TRUE)
plot_estimates <- function(data, estimates = metadata(data)[["estimates_dynamics"]],
                           delta_t = TRUE, dynamics = TRUE, distance_conditions = TRUE) {
  # bind variables to function
  r <- NULL

  # Input checks
  if (!is.list(estimates) & !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a SummarizedExperiment object or provide estimates")
  }
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    estimates <- metadata(data)[["estimates_dynamics"]]
  }
  if (!is.list(estimates)) {
    stop("'estimates' must be a list of dataframes obtained by estimates_dynamics()")
  }
  if (!is.logical(delta_t)) {
    stop("'delta_t' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(dynamics)) {
    stop("'dynamics' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(distance_conditions)) {
    stop("'distance_conditions' must be either 'TRUE' or 'FALSE'")
  }

  # binding of global variables
  `97.5%` <- NULL
  `2.5%` <- NULL
  condition <- NULL
  timepoint_1 <- NULL
  timepoint_2 <- NULL
  times <- NULL
  condition_1 <- NULL
  condition_2 <- NULL
  conditions <- NULL
  metabolite <- NULL
  time <- NULL

  plots <- list()
  plots_delta_t <- list()
  # differences between timepoints
  if (delta_t == TRUE) { # plot requested?
    if (is.data.frame(estimates[["delta_mu"]])) { # do we have a result?
      delta_mu <- estimates[["delta_mu"]]
      delta_mu <- delta_mu %>% mutate(col = ifelse(`97.5%` < 0, "CrI<0",
        ifelse(`2.5%` > 0, "CrI>0", "0inCrI")
      ))
      for (j in unique(estimates[["delta_mu"]]$condition)) {
        temp <- delta_mu %>% filter(condition == j)

        # create rank
        temp <- temp %>%
          group_by(condition, timepoint_1, timepoint_2) %>%
          arrange(mean) %>%
          mutate(r = row_number(), times = paste0(timepoint_2, "-", timepoint_1)) # delta_t = t2-t1

        for (i in unique(temp$times)) {
          temp_plot <- temp %>% filter(times == i)

          plots_delta_t[[paste0(j, "_", i)]] <-
            ggplot(temp_plot, aes(x = r, y = mean, col = col)) +
            geom_point() +
            geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2) +
            ylab("delta") +
            scale_color_manual(
              values = c("0inCrI" = "black", "CrI>0" = "green", "CrI<0" = "red")
            ) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            theme_bw() +
            xlab("metabolite") +
            scale_x_continuous(breaks = temp_plot$r, labels = temp_plot$metabolite) +
            facet_grid(cols = vars(times), rows = vars(condition)) +
            theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
            ggtitle(
              "differences between timepoints",
              "point = mean, errorbar = 95% highest density interval (CrI)"
            )
          plots[["delta_t"]] <- plots_delta_t
        }
      }
    }
    if (!is.data.frame(estimates[["delta_mu"]])) {
      message("Differences between time points can only be plotted if number of timepoints >1.")
    }
  }

  plots_distances <- list()
  # differences between timepoints
  if (distance_conditions == TRUE) { # plot requested?
    if (is.data.frame(estimates[["euclidean_distances"]])) { # do we have a result?
      distances <- estimates[["euclidean_distances"]]

      # create rank
      distances <- distances %>%
        group_by(condition_1, condition_2) %>%
        arrange(mean) %>%
        mutate(r = row_number(), conditions = paste0(condition_1, "_", condition_2)) # delta_t = t2-t1

      for (i in unique(distances$conditions)) {
        temp_plot <- distances %>% filter(conditions == i)

        plots_distances[[i]] <-
          ggplot(temp_plot, aes(y = r, x = mean)) +
          geom_point() +
          geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.2) +
          xlab("euclidean distancs between dynamics vectors") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          theme_bw() +
          ylab("metabolite") +
          scale_y_continuous(breaks = temp_plot$r, labels = temp_plot$metabolite) +
          facet_grid(cols = vars(conditions)) +
          ggtitle(
            "differences between dynamics of experimental conditons",
            "point = mean, errorbar = 95% highest density interval (CrI)"
          )
      }
      plots[["distance_conditions"]] <- plots_distances
    }
    if (!is.data.frame(estimates[["euclidean_distances"]])) {
      message("Differences of dynamics between conditons can only be plotted if number of conditions >1.")
    }
  }


  # dynamics
  if (dynamics == TRUE) {
    temp_d <- estimates[["mu"]]
    temp_d <- temp_d %>% select(metabolite, condition, time, mean)

    plots[["dynamcis"]] <-
      ggplot(temp_d, aes(
        x = as.factor(time),
        y = mean, group = metabolite, col = metabolite
      )) +
      geom_line() +
      xlab("time point") +
      scale_color_viridis_d() +
      ylab("estimated deviation from mean abundance") +
      theme_bw() +
      theme(legend.position = "none") +
      facet_grid(rows = vars(condition)) +
      ggtitle("dynamics", "color = metabolite, row label = condition")
  }

  return(plots)
}

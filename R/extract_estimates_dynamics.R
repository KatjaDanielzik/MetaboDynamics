#' extract_estimates_dynamics
#'
#' Extracts the mean concentrations (mu) at every timepoint from the dynamics model fit, the 95% highest density interval (HDI), the estimated standard deviation of metabolite concentrations at every time point (sigma), and the pooled standard deviation of every metabolite over all timepoints (lambda).
#' Additionally samples from the posterior of mu can be drawn. This can be helpful if p.e. one wants to estimate the clustering precision. Lambda can be used for clustering algorithms such as VSClust that also take the variance into account.
#'
#' @param data dataframe used for modeling
#' @param M number of metabolites, default requires a column in data named "metabolite"
#' @param t number of unique timepoints in data, the default requires a column named "time"
#' @param kegg column in "data" that contains the KEGG IDs or other identifier of metabolites
#' @param condition name of column in dataframe data that specifies the experimental condition, default requires a column named "dose"
#' @param fits list of model fits for which estimates should be extracted
#' @param iter how many iterations were used to fit the dynamics model
#' @param warmup how many warm-up iterations were used to fit the dynamics model
#' @param chains how many chains were used to fit the dynamics model
#' @param samples how many posterior samples should be drawn (p.e. for check of clustering precision)
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats runif
#'
#' @return a list of dataframes (one per experimental condition) that contains
#' the estimates at the timepoints and samples from the posterior
#' (number as specified in samples)
#' @export
#'
#' @examples
#' data("data_sim")
#' data <- data_sim[data_sim$condition == "A" & data_sim$metabolite == "ATP", ]
#' fits <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", time = "time",
#'   condition = "condition", max_treedepth = 14,
#'   adapt_delta = 0.999, iter = 4000, cores = 1, chains = 1
#' )
#' estimates <- extract_estimates_dynamics(
#'   data = data, fits = fits, iter = 4000,
#'   chains = 1, condition = "condition"
#' )
#' head(estimates)
#'
extract_estimates_dynamics <- function(data, M = length(unique(data$metabolite)),
                                       t = length(unique(data$time)),
                                       kegg = "KEGG" ,
                                       condition = "dose", fits, iter = 2000,
                                       warmup = iter / 4, chains = 4, samples = 1) {
  # bind variables
  dynamics_loc_cpc <- NULL
  temp_t <- NULL
  temp <- NULL
  metabolite.ID <- NULL
  metabolite <- NULL
  higher <- NULL
  lower <- NULL
  timepoints <- NULL

  conditions <- unique(data[[condition]])

  dynamics <- list()

  # get estimates
  for (i in names(fits))
  {
    # generate dataframe for storage
    dynamics_log_cpc <- as.data.frame(cbind(condition = i, metabolite.ID = 1:M))
    # select only models for which delta mu with 0Gy/0h was computed
    fit <- fits[[i]]
    # generate n=samples random draw numbers from posterior, but same draw
    # from every delta alpha as they are dependent on each other
    x <- floor(runif(samples, min = 1, max = (iter - warmup) * chains))

    # posterior summary
    pS <- as.data.frame(rstan::summary(fit)$summary)

    # extract posterior samples
    mu_posterior <- as.matrix(fit)


    for (m in 1:M) {
      # single loops to get reasonable order for clustering
      for (j in 1:t) {
        # get means of estimated means at timepoints
        dynamics_log_cpc[m, paste0("mu", j, "_mean")] <-
          pS[paste0("mu[", m, ",", j, "]"), "mean"]
      }
      for (j in 1:t) {
        # lower border of 95% CrI
        dynamics_log_cpc[m, paste0("mu", j, "_lower")] <-
          pS[paste0("mu[", m, ",", j, "]"), "2.5%"]
      }
      for (j in 1:t) {
        # higher border of 95% CrI
        dynamics_log_cpc[m, paste0("mu", j, "_higher")] <-
          pS[paste0("mu[", m, ",", j, "]"), "97.5%"]
      }
      for (j in 1:t) {
        # get means of sigma
        dynamics_log_cpc[m, paste0("sigma", j, "_mean")] <-
          pS[paste0("sigma[", m, ",", j, "]"), "mean"]
      }
      for (j in 1:t) {
        # lower border of 95% CrI
        dynamics_log_cpc[m, paste0("sigma", j, "_lower")] <-
          pS[paste0("sigma[", m, ",", j, "]"), "2.5%"]
      }
      for (j in 1:t) {
        # higher border of 95% CrI
        dynamics_log_cpc[m, paste0("sigma", j, "_higher")] <-
          pS[paste0("sigma[", m, ",", j, "]"), "97.5%"]
      }
      # get means of lambda
      dynamics_log_cpc[m, "lambda_mean"] <- pS[paste0("lambda[", m, "]"), "mean"]
      # lower border of 95% CrI
      dynamics_log_cpc[m, "lambda_lower"] <- pS[paste0("lambda[", m, "]"), "2.5%"]
      # higher border of 95% CrI
      dynamics_log_cpc[m, "lambda_higher"] <- pS[paste0("lambda[", m, "]"), "97.5%"]

      # delta mus
      # we have one less delta_mu than timepoints
      for (j in 1:(t - 1))
      {
        # mean
        dynamics_log_cpc[m, paste0("delta", j, j + 1, "_mean")] <-
          pS[paste0("delta_mu[", m, ",", j, "]"), "mean"]
        # lower border of 95% CrI
        dynamics_log_cpc[m, paste0("delta", j, j + 1, "_lower")] <-
          pS[paste0("delta_mu[", m, ",", j, "]"), "2.5%"]
        # higher border of 95% CrI
        dynamics_log_cpc[m, paste0("delta", j, j + 1, "_higher")] <-
          pS[paste0("delta_mu[", m, ",", j, "]"), "97.5%"]
      }


      for (k in 1:length(x))
      {
        for (j in 1:t) {
          dynamics_log_cpc[m, paste0("mu", j, "_sample", k)] <-
            mu_posterior[, paste0("mu[", m, ",", j, "]")][x[k]]
        }
      }
    }

    # get correct assignment of metabolites and metabolite.ID used in modelling
    metabolites <- as.data.frame(cbind(
      metabolite = unique(data$metabolite),
      metabolite.ID = as.numeric(as.factor(unique(data$metabolite))),
      KEGG = unique(data[[kegg]])
    ))

    dynamics_log_cpc <- left_join(dynamics_log_cpc, metabolites,
      by = "metabolite.ID"
    )
    dynamics_log_cpc <- dynamics_log_cpc %>% relocate("metabolite", .after = "metabolite.ID")
    dynamics_log_cpc <- unique(dynamics_log_cpc)
    dynamics[[i]] <- dynamics_log_cpc
  }

  # visualize
  temp <- dynamics[[1]]
  # bind if multiple conditions are analyzed
  if (length(names(dynamics)) > 1) {
    for (i in 2:length(names(dynamics))) {
      temp <- rbind(temp, dynamics[[i]])
    }
    rm(i)
  }

  # differences between timepoints
  temp_t <- temp[, c(1:3, (6 * t + 7):(6 * t + 15))]
  temp_t <- temp_t %>% pivot_longer(
    cols = -c(condition, metabolite.ID, metabolite),
    names_to = c("timepoints", ".value"), names_sep = "_"
  )
  temp_t <- temp_t %>% mutate(col = ifelse(higher < 0, "HDI>0",
    ifelse(lower > 0, "HDI<0", "0inHDI")
  ))
  dynamics[["plot_timepoint_differences"]] <-
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

  # dynamics
  temp_d <- temp[, c(1:3, 4:(t + 3))]
  temp_d <- temp_d %>% pivot_longer(
    cols = -c(condition, metabolite.ID, metabolite),
    names_to = c("timepoints", ".value"), names_sep = "_"
  )
  dynamics[["plot_dynamics"]] <-
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
  rm(temp, temp_t, temp_d)

  return(dynamics)
  rm(pS, fit, i, mu_posterior, j, k, x, conditions, metabolites, dynamics_loc_cpc)
}

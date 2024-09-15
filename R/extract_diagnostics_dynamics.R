#' Extracts diagnostic criteria from numeric fit of Bayesian model of dynamics
#'
#' gathers number of divergences, rhat values, number of effective
#' samples (n_eff) and provides plots for diagnostics criteria as well as
#' posterior predictive checks. Dataframe "model_diagnostics" contains
#' information about experimental condition, number of divergent transitions
#' and rhat and neff values for all timepoints.
#' @param data dataframe used to fit dynamics model
#' @param N number of rows in that dataframe
#' @param M number of metabolites in experimental dataset
#' @param t number of timepoints in experimental dataset
#' @param fits list of models for which diagnostics should be extracted, is the
#' object that gets returned by fit_dynamics_model()
#' @param warmup number of warmup iterations used for model fit
#' @param iter number of iterations used for model fit
#' @param chains number of chains used for model fit
#' @param scaled_measurement concentration values used to model fit, should be normalized by
#' experimental condition and metabolite to mean of zero and standard deviation
#' of one
#' 
#' @seealso [fit_dynamics_model()]/[extract_estimates_dynamics()]
#'
#' @import tidyr
#' @import ggplot2
#' @import dplyr
#' @importFrom SummarizedExperiment colData
#'
#' @return a list which contains diagnostics criteria of all conditions in a
#' dataframe (named "model_diagnostics") and one dataframe per condition that
#' contains necessary information for Posterior predictive check
#' (named "PPC_condition"). Additionally plots for diagnostics and PPC named
#' "plot_criteria" and "plot_PPC_condition" are returned.
#' @export
#'
#' @examples
#' data("data_sim")
#' # only run after fit_dynamics_model(intra): see Vignette and documentation
#' # of function
#' data_sim <- as.data.frame(SummarizedExperiment::colData(data_sim))
#' data <- data_sim[data_sim$condition == "A" & data_sim$metabolite == "ATP", ]
#' fits <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", time = "time",
#'   condition = "condition", max_treedepth = 14,
#'   adapt_delta = 0.999, iter = 4000, cores = 1, chains = 1
#' )
#' diagnostics <- extract_diagnostics_dynamics(
#'   data = data, iter = 4000, fits = fits,
#'   chains = 1, scaled_measurement = "m_scaled"
#' )
#' diagnostics[["plot_neff"]]
#' diagnostics[["plot_rhat"]]

extract_diagnostics_dynamics <- function(data, N = nrow(data),
                                         M = length(unique(data$metabolite)),
                                         t = length(unique(data$time)),
                                         iter = 2000, warmup = iter / 4, chains = 4,
                                         fits, scaled_measurement = "m_scaled") {
  
  # check input class and convert SummarizedExperiment to dataframe
  if(is(data,"SummarizedExperiment")){
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }
  
  # create list to store all subsequent results
  list_diagnostics <- list()

  # diagnostic criteria from model
  diagnostics_dynamics <- as.data.frame(matrix(ncol = 4 + 2 * t))
  names_d <- c()
  for (i in 1:t) {
    names_d[i] <- paste0("rhat_mu", i, "_mean")
  }
  for (i in (t + 1):(2 * t)) {
    names_d[i] <- paste0("neff_mu", i - t, "_mean")
  }
  colnames(diagnostics_dynamics) <- c(
    "metabolite.ID", "condition",
    "divergences", "treedepth_error", names_d
  )
  rm(i)
  # bind columns created in function to data frame
  log_cpc_stand <- NULL
  neff.mu <- NULL
  neff.value <- NULL
  rhat.mu <- NULL
  rhat.value <- NULL
  time.ID <- NULL
  treedepth_error <- NULL


  for (i in names(fits)) {
    # load fit
    fit <- fits[[i]]
    # gather diagnostic criteria
    divergences <- rstan::get_num_divergent(fit)
    max_treedepth <- rstan::get_num_max_treedepth(fit)
    rhat <- rstan::summary(fit)$summary[, "Rhat"]
    n_eff <- rstan::summary(fit)$summary[, "n_eff"]
    condition <- i

    # gather in one dataframe
    for (m in 1:M) {
      rhat_mu_mean <- as.data.frame(matrix(ncol = t))
      neff_mu_mean <- as.data.frame(matrix(ncol = t))
      start <- (m - 1) * t + 1
      end <- (m * t)
      rhat_mu_mean[, 1:t] <- rhat[start:end]
      neff_mu_mean[, 1:t] <- n_eff[start:end]
      colnames(rhat_mu_mean) <- names_d[1:t]
      colnames(neff_mu_mean) <- names_d[(t + 1):(2 * t)]
      temp <- cbind(
        metabolite.ID = m, condition = condition,
        divergences = divergences, treedepth_error = max_treedepth,
        rhat_mu_mean, neff_mu_mean
      )
      diagnostics_dynamics <- rbind(diagnostics_dynamics, temp)
    }
    rm(m)

    # Posterior predictive check
    # turn y_rep from model fit into long format and add metabolite.ID and time.ID
    posterior <- as.data.frame(fit, pars = "y_rep")
    posterior <- pivot_longer(posterior,
      names_to = "parameter",
      values_to = "posterior", cols = 1:all_of(M * t)
    )
    posterior$metabolite.ID <- as.numeric(rep(rep(1:M, t), (iter - warmup) * chains))
    posterior$time.ID <- as.factor(rep(rep(1:t, each = M), (iter - warmup) * chains))

    # assign metabolite and time id to data
    PPC <- data
    PPC$metabolite.ID <- as.numeric(as.factor(PPC$metabolite))
    PPC$time.ID <- as.factor(as.numeric(as.factor(as.numeric(PPC$time))))


    list_diagnostics[[paste0("posterior", i)]] <- posterior

    # visualize PPC
    list_diagnostics[[paste0("plot_PCC_", i)]] <-
      ggplot(posterior, aes(x = time.ID)) +
      geom_violin(aes(y = posterior), scale = "count") +
      geom_jitter(data = PPC, aes_string(x = "time.ID", y = scaled_measurement), width = 0.05) + # aes_string allows us to use predefined variables
      theme_bw() +
      ylim(-5, 5) + # we standardized data so we are not expecting much smaller or bigger values
      ggtitle(
        paste0(
          "Posterior predicitve check ", i,
          ": points within violins?"
        ),
        "violins=posterior, points=data"
      )
  }
  diagnostics_dynamics <- diagnostics_dynamics[-1, ]
  list_diagnostics[["model_diagnostics"]] <- diagnostics_dynamics
  # cleanup
  rm(i, fit, condition, n_eff, max_treedepth, divergences, rhat, PPC)

  # visualize
  list_diagnostics[["plot_divergences"]] <-
    ggplot2::ggplot(diagnostics_dynamics, aes(
      x = condition,
      y = as.numeric(divergences)
    )) +
    geom_violin() +
    xlab("condition") +
    theme_bw() +
    ylab("number of divergent transitions") +
    ggtitle("diagnostics_dynamics", "divergent transitions?")

  list_diagnostics[["plot_treedepth_error"]] <-
    ggplot2::ggplot(diagnostics_dynamics, aes(
      x = condition,
      y = as.numeric(treedepth_error)
    )) +
    geom_violin() +
    xlab("condition") +
    theme_bw() +
    ylab("number of exceeded treedepth") +
    ggtitle("diagnostics_dynamics", "maximum treedepth exceeded?")


  temp <- diagnostics_dynamics %>% pivot_longer(
    cols = 5:(5 - 1 + t),
    names_to = "rhat.mu",
    values_to = "rhat.value"
  )
  list_diagnostics[["plot_rhat"]] <-
    ggplot2::ggplot(temp, aes(
      x = as.factor(as.numeric(as.factor(rhat.mu))),
      y = as.numeric(rhat.value)
    )) +
    geom_violin() +
    theme_bw() +
    facet_wrap(~condition) +
    ylab("Rhat") +
    xlab("timepoint") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_hline(yintercept = 1.01, col = "red") +
    ggtitle("diagnostics_dynamics", "Rhat < 1.01 ?")

  # number of effective samples
  temp <- diagnostics_dynamics %>% pivot_longer(
    cols = 9:(9 - 1 + t),
    names_to = "neff.mu",
    values_to = "neff.value"
  )
  list_diagnostics[["plot_neff"]] <-
    ggplot2::ggplot(temp, aes(
      x = as.factor(as.numeric(as.factor(neff.mu))),
      y = as.numeric(neff.value)
    )) +
    geom_violin() +
    theme_bw() +
    facet_wrap(~condition) +
    ylab("number of effective samples") +
    xlab("timepoint") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_hline(yintercept = 100, col = "red") +
    ggtitle("diagnostics_dynamics", "number of effective samples >100 ?")

  return(list_diagnostics)

  # cleanup
  rm(temp)
}

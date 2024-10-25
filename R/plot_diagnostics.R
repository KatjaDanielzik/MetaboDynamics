#' Plot diagnostic criteria of numerical fit of Bayesian model of dynamics
#'
#' @param diagnostics dataframe containing diagnostics criteria from the numerical
#' fit of Bayesian model of dynamics obtained by function diagnostics_dynamics()
#' @param data dataframe or colData of a SummarizedExperiment used to fit dynamics model
#' @param divergences should number of divergent transitions be visualized?
#' @param max_treedepth should number of exeeded maximum treedepth be visualized?
#' @param Rhat should Rhat be visualized?
#' @param n_eff should number of effective samples be visualized?
#'
#' @return plots of diagnostic criteria of numerical fit of Bayesian model of 
#' dynamics
#' @seealso  parent function [diagnostics_dynamics()]
#' visualization function for posterior predictive check /[plot_PPC()]
#' @export
#'
#' @import ggplot2
#' @examples
#' data("longitudinalMetabolomics")
#' # only run after fit_dynamics_model(intra): see Vignette and documentation
#' # of function
#' longitudinalMetabolomics <- as.data.frame(SummarizedExperiment::colData(longitudinalMetabolomics))
#' data <- longitudinalMetabolomics[longitudinalMetabolomics$condition == "A" 
#'         & longitudinalMetabolomics$metabolite == "ATP", ]
#' fits <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", time = "time",
#'   condition = "condition", max_treedepth = 14,
#'   adapt_delta = 0.999, iter = 4000, cores = 1, chains = 1
#' )
#' diagnostics <- diagnostics_dynamics(
#'   data = data, iter = 4000, fits = fits,
#'   chains = 1, scaled_measurement = "m_scaled"
#' )
#' plot_diagnostics(diagnostics = diagnostics[["model_diagnostics"]], data = data)
 
plot_diagnostics <- function(diagnostics, data, divergences=TRUE, 
                             max_treedepth=TRUE,
                             Rhat=TRUE, n_eff=TRUE){
  # bind variables to function
  neff.mu <- NULL
  neff.value <- NULL
  rhat.mu <- NULL
  rhat.value <- NULL
  condition <- NULL
  treedepth_error <- NULL
  
  # input
  t <- length(unique(data$time))
  # vector for storage  
  
  plots <- list()
  # number of divergent transitions
  if(divergences==TRUE){
  plots[["divergences"]] <- ggplot2::ggplot(diagnostics, aes(
    x = condition,
    y = as.numeric(divergences)
  )) +
  geom_violin() +
  xlab("condition") +
  theme_bw() +
  ylab("number of divergent transitions") +
  ggtitle("diagnostics_dynamics", "divergent transitions?")
  }
  
  if(max_treedepth==TRUE){
  # number of interations exceeding maximum treedepth
  plots[["max_treedepth"]] <-
  ggplot2::ggplot(diagnostics, aes(
    x = condition,
    y = as.numeric(treedepth_error)
  )) +
  geom_violin() +
  xlab("condition") +
  theme_bw() +
  ylab("number of exceeded treedepth") +
  ggtitle("diagnostics_dynamics", "maximum treedepth exceeded?")
  }
  
  
  # Rhat
  if(Rhat==TRUE){
temp <- diagnostics %>% pivot_longer(
  cols = 5:(5 - 1 + t),
  names_to = "rhat.mu",
  values_to = "rhat.value"
  )
  plots[["Rhat"]] <- 
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
  }
  
  
# number of effective samples
  if(n_eff==TRUE){
temp <- diagnostics %>% pivot_longer(
  cols = 9:(9 - 1 + t),
  names_to = "neff.mu",
  values_to = "neff.value"
)
  plots[["n_eff"]]<-
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
}
}


#' Plots posterior predictive check of numerical fit of Bayesian dynamics model
#'
#' @param posterior one dataframe per condition that
#' contains necessary information for Posterior predictive check
#' obtained by function diagnostics_dynamics()(named "PPC_condition") 
#' @param data dataframe or colData of a SummarizedExperiment used to fit dynamics model
#' @param scaled_measurement column name of concentration values used to model fit, should be normalized by
#' experimental condition and metabolite to mean of zero and standard deviation
#' of one
#' @return visual posterior predictive check
#' @export
#' @seealso  parent function [diagnostics_dynamics()]
#' visualization function for diagnostics [plot_diagnostics()]
#'
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
#' plot_PPC(posterior = diagnostics[["posterior_A"]], data = data, 
#' scaled_measurement = "m_scaled")
 
 plot_PPC <- function(posterior,data,scaled_measurement){
  # bind variables to function
   time.ID <- NULL
   
  # prepare data for PPC
  # assign metabolite and time id to data
  PPC <- data
  PPC$metabolite.ID <- as.numeric(as.factor(PPC$metabolite))
  PPC$time.ID <- as.factor(as.numeric(as.factor(as.numeric(PPC$time))))
  
  plot <- ggplot(posterior, aes(x = time.ID)) +
    geom_violin(aes(y = posterior), scale = "count") +
    geom_jitter(data = PPC, aes_string(x = "time.ID", y = scaled_measurement), width = 0.05) + # aes_string allows us to use predefined variables
    theme_bw() +
    ylim(-5, 5) + # we standardized data so we are not expecting much smaller or bigger values
    ggtitle(
      paste0(
        "Posterior predicitve check ",
        ": points within violins?"
      ),
      "violins=posterior, points=data"
    )
  return(plot)
}
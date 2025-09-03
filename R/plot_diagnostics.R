#' Plot diagnostic criteria of numerical fit of Bayesian model of dynamics
#'
#' @param diagnostics dataframe containing diagnostics criteria from the numerical
#' fit of Bayesian model of dynamics obtained by function diagnostics_dynamics()
#' @param data dataframe or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used to fit dynamics model
#' must contain column "time"
#' @param assay of the \link[SummarizedExperiment]{SummarizedExperiment} object that was used to fit the dynamics
#' model
#' @param divergences should number of divergent transitions be visualized?
#' @param max_treedepth should number of exeeded maximum treedepth be visualized?
#' @param Rhat should Rhat be visualized?
#' @param n_eff should number of effective samples be visualized?
#'
#' @return plots of diagnostic criteria of numerical fit of Bayesian model of
#' dynamics
#' @seealso  parent function [diagnostics_dynamics()]
#' visualization function for posterior predictive check [plot_PPC()]
#' @export
#'
#' @import ggplot2
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite %in% c("ATP", "ADP")]
#' data <- fit_dynamics_model(
#'   model = "scaled_log",
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- diagnostics_dynamics(
#'   data = data, assay = "scaled_log",
#'   iter = 2000, chains = 1,
#'   fit = metadata(data)[["dynamic_fit"]]
#' )
#' plot_diagnostics(data = data, assay = "scaled_log")
#' 
plot_diagnostics <- function(
    data, assay = "scaled_log",
    diagnostics = metadata(data)[["diagnostics_dynamics"]][["model_diagnostics"]],
    divergences = TRUE,
    max_treedepth = TRUE,
    Rhat = TRUE, n_eff = TRUE) {
  # bind variables to function
  neff.mu <- NULL
  neff.value <- NULL
  rhat.mu <- NULL
  rhat.value <- NULL
  condition <- NULL
  treedepth_error <- NULL

  # input checks
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
    diagnostics <- metadata(data)[["diagnostics_dynamics"]][["model_diagnostics"]]
  }

  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!is.data.frame(diagnostics)) {
    stop("'diagnostics' must be a dataframe obtained by diagnostics_dynamics()")
  }
  if (!is.logical(divergences)) {
    stop("'divergences' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(max_treedepth)) {
    stop("'max_treedepth' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(Rhat)) {
    stop("'Rhat' must be either 'TRUE' or 'FALSE'")
  }
  if (!is.logical(n_eff)) {
    stop("'n_eff' must be either 'TRUE' or 'FALSE'")
  }
  if (!all(c("time") %in% colnames(data_df))) {
    stop("'data' must contain a column named 'time'")
  }


  # input
  t <- length(unique(data_df$time))
  # vector for storage

  plots <- list()
  # number of divergent transitions
  if (divergences == TRUE) {
    plots[["divergences"]] <- ggplot2::ggplot(diagnostics, aes(
      x = condition,
      y = divergences
    )) +
      geom_violin() +
      xlab("condition") +
      theme_bw() +
      ylab("number of divergent transitions") +
      ggtitle("diagnostics_dynamics", "divergent transitions?")
  }

  if (max_treedepth == TRUE) {
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
  if (Rhat == TRUE) {
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
  if (n_eff == TRUE) {
    temp <- diagnostics %>% pivot_longer(
      cols = (5 + t):(5 - 1 + 2 * t),
      names_to = "neff.mu",
      values_to = "neff.value"
    )
    plots[["n_eff"]] <-
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
  return(plots)
}

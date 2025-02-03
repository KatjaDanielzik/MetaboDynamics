#' Plot results of over-representation analysis with ORA_hypergeometric()
#'
#' @param data result dataframe from ORA_hypergeometric() or \link{SummarizedExperiment}
#' object where the ORA_hypergeometric() results are stored in metadata(data)
#' under "ORA_tested_column"
#' @param tested_column KEGG module hierarchy level on which ORA was executed
#' @return a plot of the over-representation analysis
#' @export
#' @import ggplot2
#' @import dplyr
#' @seealso do over-represenation analysis of KEGG functional modules [ORA_hypergeometric()]
#' @examples
#' data("longitudinalMetabolomics")
#' data("modules_compounds")
#' head(modules_compounds)
#' data("metabolite_modules")
#' head(metabolite_modules)
#' # middly hierachy
#' longitudinalMetabolomics <- ORA_hypergeometric(
#'   data = longitudinalMetabolomics,
#'   annotations = metabolite_modules,
#'   background = modules_compounds,
#'   tested_column = "middle_hierarchy"
#' )
#' plot_ORA(longitudinalMetabolomics)
plot_ORA <- function(data, tested_column = "middle_hierarchy") {
  # bind variables to function
  OvE_gen <- NULL
  module_name <- NULL
  OvE_gen_lower <- NULL
  OvE_gen_higher <- NULL
  OvE_gen_median <- NULL
  cluster <- NULL
  condition <- NULL

  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  if (is(data, "SummarizedExperiment")) {
    a_clusters <- metadata(data)[[paste0("ORA_", tested_column)]]
  }
  if (is(data, "data.frame")) {
    a_clusters <- data
  }

  # get module name (in ORA_hypergeometric this is the 3rd column of the result)
  module_name <- colnames(a_clusters)[3]
  module_name <- as.symbol(module_name)
  module_name <- enquo(module_name)

  # reduce rows for faster plotting
  a_clusters <- a_clusters %>% select(
    condition, cluster, !!module_name, OvE_gen,
    OvE_gen_median, OvE_gen_lower,
    OvE_gen_higher
  )
  a_clusters <- unique(a_clusters)

  # color code for visualization
  a_clusters <- a_clusters %>% mutate(col = ifelse(log(OvE_gen_higher) < 0,
    "ICR<0",
    ifelse(log(OvE_gen_lower) > 0,
      "ICR>0", "ICR includes 0"
    )
  ))

  plot <- ggplot(
    a_clusters,
    aes(x = log(as.numeric(OvE_gen)), y = (!!module_name), col = col)
  ) +
    geom_errorbarh(aes(
      xmin = log(as.numeric(OvE_gen_lower)),
      xmax = log(as.numeric(OvE_gen_higher))
    )) +
    geom_point(aes(x = log(as.numeric(OvE_gen_median)))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw() +
    scale_color_manual(values = c("ICR includes 0" = "black", "ICR>0" = "green", "ICR<0" = "red")) +
    xlab("log(p(OvE))") +
    guides(col = "none") +
    facet_grid(cols = vars(cluster), rows = vars(condition)) +
    ggtitle(
      "hypergeometric ORA",
      "median and 95% interquantile range, panels=clusterID"
    )
  return(plot)
}

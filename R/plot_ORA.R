#' Plot results of over-representation analysis with ORA_hypergeometric()
#'
#' @param ORA result dataframe from 
#'
#' @return a plot of the over-representation analysis 
#' @export
#' @import ggplot2
#' @seealso [ORA_hypergeometric()]
#' @examples
#'  data("cluster")
#' head(cluster)
#' data("modules_compounds")
#' head(modules_compounds)
#' data("metabolite_modules")
#' head(metabolite_modules)
#' ORA <- ORA_hypergeometric(
#'   background = modules_compounds,
#'   annotations = metabolite_modules, clusters = cluster,
#'   tested_column = "middle_hierarchy"
#' )
#' plot_ORA(ORA)

plot_ORA <- function(ORA){
plot <- ggplot(ORA,
  aes(x = log(OvE_gen), y = module_name, col = col)
) +
  geom_errorbarh(aes(xmin = log(OvE_gen_lower), xmax = log(OvE_gen_higher))) +
  geom_point(aes(x = log(OvE_gen_median))) +
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

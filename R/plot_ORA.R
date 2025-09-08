#' Plot results of over-representation analysis with ORA_hypergeometric()
#'
#' @param data result dataframe from [ORA_hypergeometric()] or \link[SummarizedExperiment]{SummarizedExperiment}
#' object where the ORA_hypergeometric() results are stored in metadata(data)
#' under "ORA_tested_column"
#' @param tested_column KEGG module hierarchy level on which ORA was executed
#' @param patchwork should result be patchworked to results of [plot_cluster()]?
#' @param plot_cluster if patchwork = TRUE this needs to be the result of [plot_cluster()]
#' @return a plot of the over-representation analysis and lsit of plots suitable to
#' patchwork with cluster visualization if patchwork=TRUE
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
#' data("IDs")
#' head(IDs)
#' # middly hierachy
#' longitudinalMetabolomics <- ORA_hypergeometric(
#'   data = longitudinalMetabolomics,
#'   annotations = metabolite_modules,
#'   background = modules_compounds,
#'   tested_column = "middle_hierarchy",
#'   IDs = IDs
#' )
#' plot_ORA(longitudinalMetabolomics)
plot_ORA <- function(data, tested_column = "middle_hierarchy",
                     patchwork = FALSE,
                     plot_cluster = NULL) {
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

  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    a_clusters <- data
  }
  if (!is.logical(patchwork)) {
    stop("'patchwork' must be either 'TRUE' or 'FALSE'")
  }
  if (isTRUE(patchwork)) {
    if(!inherits(plot_cluster,"list"))
    stop("if 'patchwork is TRUE, plot_cluster must be provided")
  }

  # get module name (in ORA_hypergeometric this is the 3rd column of the result)
  module_name <- colnames(a_clusters)[3]
  module_name <- as.symbol(module_name)
  module_name <- enquo(module_name)
  tested_column <- as.symbol(tested_column)
  tested_column <- enquo(tested_column)

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
    scale_color_manual(
      values = c("black", "green", "red"),
      labels = c("0 in ICR", "ICR>0", "ICR<0"), name = ""
    ) +
    xlab("log(p(OvE))") +
    facet_grid(cols = vars(cluster), rows = vars(condition)) +
    ggtitle(
      "hypergeometric ORA",
      "median and 95% interquantile range, panels=clusterID"
    )
  
  ora_patchwork <- list()
  if(patchwork == TRUE){
    for(i in unique(a_clusters$condition)){
      temp <- a_clusters%>%filter(condition==i)
      # get cluster order
      cluster_order <- plot_cluster$cluster_order[[i]]
      # order cluster in ORA hypergeometric
      temp$cluster <- as.factor(temp$cluster)
      temp$cluster <- factor(temp$cluster, levels=cluster_order)
      temp <- temp[order(temp$cluster),]
      plots <- list()
      for (j in unique(temp$cluster)){
        temp_plot <- temp%>%filter(cluster==j)
        plots[[j]] <-
          ggplot(temp_plot,aes(x=OvE_gen_median,y=cluster,col=col))+
            geom_point()+
            geom_errorbarh(aes(xmin=OvE_gen_lower,xmax=OvE_gen_higher))+
            facet_grid(cols=vars(!!tested_column))+
            theme_bw()+
            geom_vline(xintercept = 0, linetype = "dashed") +
            ylab("cluster ID")+
            xlab("") +
            scale_color_manual(
            values = c("black", "green", "red"),
            labels = c("0 in ICR", "ICR>0", "ICR<0"), name = "")
        plots[[paste0(j,"space")]] <- plot_spacer() 
      }
      heights <- plot_cluster$cluster_heights[[i]] ## needs to be sorted!
      p <- Reduce("/",plots)
      ora_patchwork[[i]] <- p + plot_layout(heights=heights)+plot_annotation("ORA of KEGG functional modules",
                                                                         "mean and 95% ICR")
      }
      }
  return(list(plot=plot,ora_patchwork=ora_patchwork))
}

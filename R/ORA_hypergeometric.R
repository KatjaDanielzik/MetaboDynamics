#' OverRepresentationAnalysis with a hypergeometric model
#'
#' Testing the hypothesis that certain KEGG modules are over-represented in
#' clusters of metabolites.
#' A module is considered over-represented in a cluster the number of
#' metabolites in a cluster being annotated to a functional module (n_obs)
#' is higher than the expected number of metabolites in a cluster of this size
#' being annotated to a functional module (n_theo).
#' We can calculate the OvE (Observed versus Expected=n_obs/n_theo) and show the
#' probabilites of these ratios.
#' log(p(OvE))>0 indicates an over-representation of the functional module in
#' the cluster, log(p(OvE))<0 an under-representation.
#'
#' @param background dataframe that contains which metabolites
#' (represented as KEGG ID) are annotated to functional modules in general
#' @param annotations to which functional modules our experimental metabolites
#' are annotated
#' @param clusters dataframe containing columns "metabolite","cluster" and a
#' column specifying the experimental condition called "condition"
#' @param tested_column column that is in background and annotations and on
#' which the hypergeometric model shall be executed
#
#' @return a list with a dataframe containing the ORA results and a plot of
#' the ORA
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats rhyper
#'
#' @examples
#' data("cluster")
#' # middly hierachy
#' ORA <- ORA_hypergeometric(background = modules_compounds, annotations = metabolite_modules, clusters = cluster, tested_column = "middle_hierarchy")
#' ORA[["plot_ORA"]]
#' # lower hierachy
#' ORA_lower <- ORA_hypergeometric(background = modules_compounds, annotations = metabolite_modules, clusters = cluster[cluster$condition=="A",], tested_column = "lower_hierarchy")
#' ORA_lower[["plot_ORA"]]


ORA_hypergeometric <- function(background, annotations,
                               clusters, tested_column = "middle_hierarchy") {
  # return object
  ORA_hyper <- list()

  # attach new variables to function
  N_b <- NULL
  N <- NULL
  M <- NULL
  a_clusters <- NULL
  mapped_m <- NULL
  a <- NULL
  temp_N <- NULL
  a_list <- NULL
  series <- NULL
  hits_in_module <- NULL
  background_module <- NULL
  total_background <- NULL
  total_in_cluster <- NULL
  n_gen <- NULL
  n_obs <- NULL
  OvE_gen <- NULL
  median <- NULL
  OvE_gen_lower <- NULL
  OvE_gen_higher <- NULL
  OvE_gen_median <- NULL
  module <- NULL
  module_name <- NULL
  condition <- NULL
  cluster <- NULL

  # all unique metabolites in Background -> uniquely to avoid bias p.e.
  # for side-compounds
  N_b <- unique(background$kegg_id)
  # all experimental metabolites that are mapped to KEGG modules
  N <- unique(annotations[!is.na(annotations$module_id), ]$KEGG)
  # list of all middle hierachy modules
  M <- unique(background[[tested_column]])
  ## get dataframe with only experimental metabolites that are mapped to
  # KEGG module
  mapped_m <- annotations[annotations$KEGG %in% N, ]

  # internal helper function to retrieve annotated KEGG IDs of a module
  #' @keywords internal
  .get_module_background <- function(M) {
    return(background[background[tested_column] == M, ]$kegg_id)
  }

  # get list of Modules and corresponding KEGG IDs
  a_b_list <- sapply(M, .get_module_background)
  # extract number of unique background metabolites in module
  a_b <- c()
  for (i in 1:length(M))
  {
    # extract number of unique experimental metabolites in module
    a_b[i] <- length(unique(unlist(a_b_list[names(a_b_list)][M[i]])))
  }
  rm(i, a_b_list)


  # retrieve cluster membership information of experimental metabolites
  a_clusters <- as.data.frame(matrix(ncol = 5))
  colnames(a_clusters) <- c(
    "condition", "cluster", "module",
    "total_in_cluster", "hits_in_module"
  )
  ## get list of experimental metabolites in modules
  #' internal helper function to retrieve modules experimental metabolites
  #' are annotated to
  #' @keywords internal
  .get_module <- function(M) {
    return(temp[temp[tested_column] == M, ]$KEGG)
  }

  for (j in unique(clusters$condition)) {
    temp1 <- left_join(mapped_m, clusters[clusters$condition == j,],
      by = "metabolite"
    )
    for (i in 1:length(unique(temp1$cluster))) {
      temp <- temp1[temp1$cluster == i, ]
      a_list <- sapply(M, .get_module)
      # vector of metabolites in this cluster sample
      N <- unique(temp$KEGG)
      ## extract number of unique metabolites
      a <- c(1:length(M))
      for (k in 1:length(M))
      {
        # extract number of unique experimental metabolites in module
        a[k] <- length(unique(unlist(a_list[names(a_list)][M[k]])))
      }
      temp_N <- as.data.frame(cbind(
        condition = j, cluster = i,
        module = 1:length(M),
        total_in_cluster = length(N),
        hits_in_module = a
      ))
      a_clusters <- rbind(a_clusters, temp_N)
    }
  }
  rm(i, j, temp, temp_N, temp1, a, k, a_list, mapped_m)
  a_clusters <- a_clusters[!is.na(a_clusters$condition), ]
  a_clusters$module <- as.numeric(a_clusters$module)

  test <- as.data.frame(cbind(
    module = as.numeric(1:length(M)),
    background_module = a_b, module_name = M
  ))
  rm(M, N, N_b)
  test$module <- as.numeric(test$module)
  a_clusters <- left_join(a_clusters, test, by = "module")
  a_clusters$total_background <- sum(a_b)
  rm(test, a_b)


  # hypergeometric test
  ## get numericals
  a_clusters$total_in_cluster <- as.numeric(a_clusters$total_in_cluster)
  a_clusters$hits_in_module <- as.numeric(a_clusters$hits_in_module)
  a_clusters$background_module <- as.numeric(a_clusters$background_module)
  a_clusters$total_background <- as.numeric(a_clusters$total_background)

  ## generate column with values from 0-total_in_cluster (N):
  temp <- unique(a_clusters[, c("condition", "cluster", "total_in_cluster")])
  series <- as.data.frame(matrix(ncol = 4))
  colnames(series) <- c("condition", "cluster", "total_in_cluster", "n_theo")

  for (i in 1:nrow(temp)) {
    temp2 <- temp[i, ]
    temp3 <- cbind(temp2, n_theo = seq(from = 1, to = 100))
    series <- rbind(series, temp3)
  }
  rm(temp, temp2, temp3, i)
  series <- series[-1, ]
  a_clusters <- left_join(a_clusters, series, by = join_by(
    "condition",
    "cluster",
    "total_in_cluster"
  ), relationship = "many-to-many")
  rm(series)

  # get probability of n_theo
  # m= white balls in urn = metabolites in modules in background
  # (background_module)
  # x= drawn white balls = theoretically drawn metabolites (n_theo)
  # n= number of fails in urn= total in background-metabolites in modules
  # in background
  # k= number of draws from the urn= number of cluster metabolites
  # (total_in_cluster)

  # n_obs=0 causes problems down the line

  # add pseudocounts if n_obs=0 to avoid very small OvE
  a_clusters <- a_clusters %>% mutate(n_obs = ifelse(hits_in_module == 0, 1,
    hits_in_module
  ))

  a_clusters <- a_clusters %>%
    group_by(condition, cluster, module) %>%
    mutate(n_gen = rhyper(100,
      m = background_module,
      n = total_background - background_module,
      k = total_in_cluster
    ))

  a_clusters <- a_clusters %>% mutate(n_gen = ifelse(n_gen == 0, 1, n_gen))

  a_clusters <- a_clusters %>% mutate(
    OvE_gen = n_obs / n_gen,
    OvE_gen_lower = quantile(OvE_gen, probs = 0.025),
    OvE_gen_higher = quantile(OvE_gen, probs = 0.975),
    OvE_gen_median = median(OvE_gen)
  )

  # color code for visualization
  a_clusters <- a_clusters %>% mutate(col = ifelse(log(OvE_gen_higher) < 0,
                                                   "ICR<0",
                                                   ifelse(log(OvE_gen_lower) > 0,
                                                   "ICR>0", "ICR includes 0")))

  ORA_hyper[["ORA"]] <- a_clusters
  ORA_hyper[["plot_ORA"]] <- ggplot(
    a_clusters,
    aes(x = log(OvE_gen), y = module_name, col = col)) +
    geom_errorbarh(aes(xmin = log(OvE_gen_lower), xmax = log(OvE_gen_higher))) +
    geom_point(aes(x = log(OvE_gen_median))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw() +
    scale_color_manual(values = c("ICR includes 0" = "black", "ICR>0" = "green","ICR<0"="red")) +
    xlab("log(p(OvE))") +
    guides(col = "none") +
    facet_grid(cols = vars(cluster), rows = vars(condition)) +
    ggtitle("hypergeometric ORA",
            "median and 95% interquantile range, panels=clusterID")
  return(ORA_hyper)
}

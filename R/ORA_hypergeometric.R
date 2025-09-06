#' OverRepresentationAnalysis with a hypergeometric model
#'
#' Testing the hypothesis that certain KEGG modules are over-represented in
#' clusters of metabolites.
#' A module is considered over-represented in a cluster the number of
#' metabolites in a cluster being annotated to a functional module (n_obs)
#' is higher than the expected number of metabolites in a cluster of this size
#' being annotated to a functional module (n_theo).
#' We can calculate the OvE (Observed versus Expected = n_obs/n_theo) and show the
#' probabilities of these ratios.
#' log(p(OvE))>0 indicates an over-representation of the functional module in
#' the cluster, log(p(OvE))<0 an under-representation.
#' @seealso 
#' Function to visualize ORA results [plot_ORA()]
#' @param background dataframe that contains
#' KEGG IDs of metabolites that are assigned to functional modules, is incorporated
#' in the package [modules_compounds]
#' @param annotations dataframe tha contains information to which functional
#' modules our experimental metabolites are annotated in KEGG, can be constructed
#' by filtering the provided KEGG background [modules_compounds] for the experimental
#' metabolites
#' @param data data frame containing columns "KEGG" specifying the KEGG
#' identifiers of metabolites, "cluster" specifying the cluster ID of metabolites and a
#' column specifying the experimental condition called "condition" or if data
#' is a SummarizedExperiment or a \link[SummarizedExperiment]{SummarizedExperiment} clustering
#' solution must be stored in metadata(data)
#' under "cluster"
#' @param tested_column column that is in background and annotations and on
#' which the hypergeometric model will be executed
#
#' @return a dataframe containing the ORA results or if data is SummarizedExperiment \link[SummarizedExperiment]{SummarizedExperiment}
#' object the output is stored in metadata(data) under "ORA_tested_column"
#' @export
#'
#' @import dplyr
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats rhyper
#'
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
#' S4Vectors::metadata(longitudinalMetabolomics)[["ORA_middle_hierarchy"]]
#' 
ORA_hypergeometric <- function(background,
                               annotations,
                               data, tested_column = "middle_hierarchy") {
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
  median <- NULL
  module <- NULL
  OvE_gen <- NULL
  module_name <- NULL
  OvE_gen_lower <- NULL
  OvE_gen_higher <- NULL
  OvE_gen_median <- NULL
  cluster <- NULL
  condition <- NULL
  KEGG <- NULL
  module_id <- NULL

  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object obtained by cluster_dynamics")
  }
  if (is(data, "SummarizedExperiment")) {
    data <- metadata(data)[["cluster"]]
    # combine data of list elements
    data <- bind_rows(lapply(data_df, function(x){
      return(x$data)
    }))
  }

  if (is(data, "list")) {
    # combine data of list elements
    data <- bind_rows(lapply(data, function(x){
      return(x$data)
    }))
  }
  # convert potential tibbles into data frame
  if (is(data_df, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }

  # input checks
  if (!is.data.frame(background)) {
    stop("'background' must be a dataframe obtained by get_ORA_annotations()")
  }
  if (!is.data.frame(annotations)) {
    stop("'annotations' must be a dataframe obtained by get_ORA_annotations()")
  }
  if (!is.character(tested_column)) {
    stop("'tested_column' must be a character vector")
  }
  if (!tested_column %in% colnames(background)) {
    stop("'tested_column' must be a column of 'background'")
  }
  if (!tested_column %in% colnames(annotations)) {
    stop("'tested_column' must be a column of 'annotations'")
  }
  if (!all(c("module_id", "module_name", "KEGG", tested_column) %in% colnames(background))) {
    stop("'background' must contains columns 'module_id', 'module_name', 'KEGG' and 'tested_column'
         suitable dataframes can be obtained with get_ORA_annoations()")
  }
  if (!all(c("module_id", "module_name", "KEGG", tested_column) %in% colnames(annotations))) {
    stop("'annotations' must contains columns 'module_id', 'module_name', 'KEGG' and 'tested_column'
         suitable dataframes can be obtained with get_ORA_annoations()")
  }
  if (!all(c("KEGG", "condition", "cluster") %in% colnames(data_df))) {
    stop("'data' must contains columns 'KEGG', 'cluster' and 'condition'")
  }

  name <- tested_column

  # all experimental metabolites that are mapped to KEGG modules
  N <- unique(annotations[!is.na(annotations$module_id), ]$KEGG)
  ## get data frame with only experimental metabolites that are mapped to
  # KEGG module
  mapped_m <- annotations[annotations$KEGG %in% N, ] %>%
    distinct(KEGG, module_id, .keep_all = TRUE)

  # retrieve cluster membership information of experimental metabolites
  a_clusters <- data_df %>%
    left_join(mapped_m, join_by("KEGG"), relationship = "many-to-many") %>%
    filter(!is.na(module_id))
  # total number of metabolites in one cluster
  total_in_cluster <- a_clusters %>%
    group_by(condition, cluster) %>%
    summarise(total_in_cluster = n_distinct(KEGG))
  a_clusters <- left_join(a_clusters, total_in_cluster, join_by(condition, cluster))
  # number of cluster members in one module (selected hierarchy level)
  # get column name from user
  tested_column <- as.symbol(tested_column)
  tested_column <- enquo(tested_column)
  hits_in_module <- a_clusters %>%
    group_by(condition, cluster, (!!tested_column)) %>%
    summarise(hits_in_module = n_distinct(KEGG))
  a_clusters <- left_join(a_clusters, hits_in_module, join_by(condition, cluster, !!tested_column))
  # number of background metabolites per module
  background_module <- background %>%
    group_by(!!tested_column) %>%
    summarise(background_module = n_distinct(KEGG))
  a_clusters <- left_join(a_clusters, background_module, join_by(!!tested_column))
  # total metabolites in background
  a_clusters$total_background <- as.numeric(n_distinct(background$KEGG))

  ## generate column with values from 0-total_in_cluster (N) = theoretical possible
  # number of hits in module:
  a_clusters <- unique(a_clusters %>% select(
    condition, cluster, !!tested_column,
    total_in_cluster, hits_in_module,
    background_module, total_background
  ))
  temp <- unique(a_clusters %>% select(condition, cluster, total_in_cluster))
  temp <- temp %>%
    group_by(condition, cluster) %>%
    slice(rep(1, total_in_cluster))
  temp <- temp %>% mutate(n_theo = seq_len(unique(total_in_cluster)))
  a_clusters <- left_join(a_clusters, temp,
    join_by(condition, cluster, total_in_cluster),
    relationship = "many-to-many"
  )

  # add pseudocounts if n_obs=0 to avoid very small OvE
  a_clusters <- a_clusters %>% mutate(n_obs = ifelse(hits_in_module == 0, 1,
    hits_in_module
  ))

  # generate random numbers from hypergeometric distribution using lapply()
  a_clusters <- a_clusters %>%
    group_by(condition, cluster, !!tested_column) %>%
    mutate(n_gen = list(rhyper(1e4,
      m = background_module,
      n = total_background - background_module,
      k = total_in_cluster
    ))) %>%
    unnest(n_gen) %>%
    mutate(n_gen = ifelse(n_gen == 0, 1, n_gen)) %>%
    group_by(condition, cluster, !!tested_column) %>%
    mutate(
      OvE_gen = n_obs / n_gen,
      OvE_gen_lower = quantile(OvE_gen, probs = 0.025, na.rm = TRUE),
      OvE_gen_higher = quantile(OvE_gen, probs = 0.975, na.rm = TRUE),
      OvE_gen_median = median(OvE_gen, na.rm = TRUE)
    )
  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    name <- paste0("ORA_", name)
    metadata(data)[[name]] <- a_clusters
    return(data)
  } else {
    # otherwise, return the list of fits
    return(a_clusters)
  }
}

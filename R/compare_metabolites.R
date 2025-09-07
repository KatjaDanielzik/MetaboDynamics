#' Comparison of metabolite sets between dynamics clusters of different experimental conditions
#'
#' Uses the Jaccard Index to compare metabolite names between dynamics clusters of
#' different experimental conditions
#'
#' @param data result of [cluster_dynamics()] function: either a list of data frames or a SummarizedExperiment object
#' @return a data frame of Jaccard indices between data or if data input
#' was a SummarizedExperiment results are stored in metadata(data) under
#' "comparison_metabolites"
#'
#' @export
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @seealso Visualization of metabolite similarity [heatmap_metabolites()]/
#' compare dynamics of clusters [compare_dynamics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' longitudinalMetabolomics <- compare_metabolites(
#'   data = longitudinalMetabolomics
#' )
#' S4Vectors::metadata(longitudinalMetabolomics)[["comparison_metabolites"]]
compare_metabolites <- function(data) {
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object obtained by 'cluster_dynamics()'")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- metadata(data)[["cluster"]]
    # combine data of list elements
    data_df <- bind_rows(lapply(data_df, function(x){
      return(x$data)
    }))
  }
  if (is(data, "list")) {
    # combine data of list elements
    data_df <- bind_rows(lapply(data, function(x){
      return(x$data)
    }))
  }
  #convert potential tibbles into data frame
  if (is(data_df, "tbl")) {
    data <- as.data.frame(data_df)
  }
  if (is(data_df, "data.frame")) {
    data_df <- data_df
  }
  if (is(data_df, "tbl")) {
    data_df <- as.data.frame(data_df)
  }
  
  if (!all(c("cluster","metabolite","condition") %in% colnames(data_df))) {
    stop("'data' must contain columns named 'cluster','metabolite' and 'condition'")
  }

  # order data_df according to condition and cluster
  data_df$cluster <- as.numeric(data_df$cluster)
  data_df <- data_df[order(data_df$condition,data_df$cluster),]
  

  # Extract unique condition-cluster combinations
  unique_combinations <- unique(data_df[, c("condition", "cluster")])
  n_combinations <- nrow(unique_combinations)

  # Generate upper triangle of pairwise combinations (i.e., no self-comparisons)
  comparison_pairs <- combn(seq_len(n_combinations), 2)

  # Assign names to distances
  comparison_names <- apply(comparison_pairs, 2, function(idx) {
    paste(
      paste(unique_combinations$condition[idx[1]], unique_combinations$cluster[idx[1]], sep = "_"),
      "vs",
      paste(unique_combinations$condition[idx[2]], unique_combinations$cluster[idx[2]], sep = "_")
    )
  })

  # Function to calculate Jaccard index
  calculate_jaccard <- function(group_a, group_b) {
    .similarity(group_a, group_b) # Assuming .similarity calculates Jaccard Index
  }

  # Prepare metabolite lists for all combinations
  metabolite_lists <- lapply(seq_len(n_combinations), function(i) {
    subset <- data_df[data_df$condition == as.character(unique_combinations[i, "condition"]) &
      data_df$cluster == as.numeric(unique_combinations[i, "cluster"]), "metabolite"]
    unique(subset)
  })

  # Vectorized computation of Jaccard indices
  jaccard_indices <- apply(comparison_pairs, 2, function(idx) {
    calculate_jaccard(metabolite_lists[[idx[1]]], metabolite_lists[[idx[2]]])
  })

  # Create result data frame
  result <- data.frame(
    comparison = comparison_names,
    cluster_a = apply(comparison_pairs, 2, function(idx) {
      paste(unique_combinations[idx[1], "condition"], unique_combinations[idx[1], "cluster"], sep = "_")
    }),
    cluster_b = apply(comparison_pairs, 2, function(idx) {
      paste(unique_combinations[idx[2], "condition"], unique_combinations[idx[2], "cluster"], sep = "_")
    }),
    Jaccard = jaccard_indices
  )
  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["comparison_metabolites"]] <- result
    return(data)
  } else {
    # otherwise, return the result
    return(result)
  }
}

# Function to calculate all pairwise distances between two groups
#' compare_dynamics()
#'
#' @param group_a dataframe of one cluster of one condition
#' @param group_b dataframe of one cluster of a differenct condition than group_a
#' @param dynamics character vector specifying the columns that hold dynamic
#' estimates in data
#' @return matrix of pairwise euclidean distances between two groups of vectors
#' @keywords internal
.calculate_distances <- function(group_a, group_b, dynamics) {
  group_a <- as.matrix(group_a[, dynamics, drop = FALSE])
  group_b <- as.matrix(group_b[, dynamics, drop = FALSE])
  outer(seq_len(nrow(group_a)), seq_len(nrow(group_b)), Vectorize(function(i, j) {
    .eu(group_a[i, ], group_b[j, ])
  }))
}

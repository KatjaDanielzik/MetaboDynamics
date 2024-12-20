#' Function to calculate Jaccard index on two character vectors of metabolite
#' names
#'
#' @param group_a group of clusters of metabolites
#' @param group_b group of clusters of metabolites
#'
#' @return the Jaccard index
#' @keywords internal

.calculate_jaccard <- function(group_a, group_b) {
  .similarity(group_a, group_b)
}

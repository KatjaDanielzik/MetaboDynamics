#' euclidean distance
#' compare_dynamics()
#' @param a a numeric vector
#' @param b a numeric vector of same length as a
#' @importFrom stats dist
#' @return euclidean distance between vectors
#' @keywords internal

eu <- function(a, b) {
  temp <- rbind(a, b)
  dist <- stats::dist(temp, method = "euclidean")
  return(dist)
}
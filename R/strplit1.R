#' Title
#'
#' @param x a delimiter
#' @param split a string
#'
#' @return A character vector.
#' @export
#'
#' @examples
strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}


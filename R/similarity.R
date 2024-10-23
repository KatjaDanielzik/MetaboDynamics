#' Jaccard Index: intersection/union
#' compare_metabolites()
#' @param a a vector
#' @param b a vector
#'
#' @return Jaccard Index of a and b
#' @keywords internal
#'
.similarity <- function(a, b) {
  # test which vector is bigger and compare smaller to bigger
  if(length(a)<length(b)){
  # intersection
  temp <- is.element(a,b)==TRUE
  }
  else{
  temp <- is.element(b,a)==TRUE
  }
  intersection <- length(temp[temp==TRUE])
  # union=unique metabolites per set + intersection
  sim <- intersection / sum((length(a) - intersection), 
                            (length(b) - intersection), intersection)
  return(sim)
}
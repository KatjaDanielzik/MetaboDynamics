#' internal helper function to retrieve experimental KEGG IDs annotated to module
#' get_ORA_annotations()
#' @param M module name
#'
#' @return KEGG IDs of experimental metabolites annotated to module
#' @keywords internal
.get_module <- function(M) {
  return(temp[temp[tested_column] == M, ]$KEGG)
}
#' internal helper function to retrieve background KEGG IDs of a module
#' get_ORA_annotations()
#' @param M module name
#'
#' @return KEGG IDs of background metabolites annotated to module
#' @keywords internal
.get_module_background <- function(M) {
  return(background[background[tested_column] == M, ]$kegg_id)
}
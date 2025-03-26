#' Retrieve background and annotation information for over-representation analysis
#' (ORA)
#'
#' Uses the package KEGGREST to retrieve background and annotation information
#' needed for over-representation analysis. As KEGGREST only allows 10 queries
#' per second this might take some time to run, depending on the size of the
#' dataset and organism. The user should check afterwards if all functional modules are
#' applicable for the analysis of the dataset (p.e. organism, tissue).
#'
#' @param data data frame to be analyzed with ORA. Must at least contain a column
#'                      with KEGG IDs or a \link[SummarizedExperiment]{SummarizedExperiment}  where the metabolite
#'                      names or IDs are stored in colData
#' @param kegg column name of "data" that holds KEGG IDs
#' @param metabolite_name column name of "data" that holds metabolite names
#' @param update_background logical. Should the background information be updated?
#' Should be set to TRUE of this is the first time using this function.
#' If TRUE this may take some time.
#'
#' @seealso Dor over-representation analysis of KEGG functional modules [ORA_hypergeometric()]
#'
#' @return a list with dataframes "background" and "annotation" needed for ORA,
#' if data is a SummarizedExperiment \link[SummarizedExperiment]{SummarizedExperiment}  object annotations are stored in metadata(data)
#' under "KEGG_annotations"
#' @export
#' @import SummarizedExperiment
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGREST keggList
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#' @importFrom stats na.omit
#' @import dplyr
#' @import stringr
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite == "ATP"]
#' data <- get_ORA_annotations(
#'   data = data,
#'   kegg = "KEGG",
#'   metabolite_name = "metabolite",
#'   update_background = FALSE
#' )
#' S4Vectors::metadata(data)[["KEGG_annotations"]]
get_ORA_annotations <- function(data, kegg = "KEGG",
                                metabolite_name = "metabolite",
                                update_background = TRUE) {
  # input checks
  if (!is.logical(update_background)) {
    stop("'update_background' must be either 'TRUE' or 'FALSE'")
  }
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  if (is(data, "SummarizedExperiment")) {
    data_df <- as.data.frame(colData(data))
  }
  
  # convert potential tibbles into data frame
  if(is(data,"tbl")){
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  if (!is.character(kegg)) {
    stop("'kegg' must be a character vector specifying a column name of data")
  }
  if (!is.character(metabolite_name)) {
    stop("'metabolite_name' must be a character vector specifying a column name of data")
  }
  if (!all(c(kegg, metabolite_name) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'kegg', and 'metabolite_name'")
  }

  message("This request may take some time, please be patient.")
  # bind variables to function
  ORA_dataframes <- list()
  metabolite_modules <- NULL
  temp2 <- NULL
  temp <- NULL
  temp3 <- NULL
  modules_compounds <- NULL

  metabolite_modules <- unique(as.data.frame(cbind(
    metabolite = data[[metabolite_name]],
    KEGG = data[[kegg]]
  )))

  # use lapply() to apply KEGGREST::keggGet() to each unique KEGG ID in parallel
  temp2 <- do.call(rbind, lapply(unique(data[[kegg]]), function(i) {
    if (!is.na(i)) {
      temp <- tryCatch(
        {
          Sys.sleep(0.2) # add a delay of 1 second between requests
          KEGGREST::keggGet(i)
        },
        error = function(e) NULL
      )
      temp <- as.data.frame(temp[[1]]$MODULE)
      temp <- as.data.frame(cbind(
        KEGG = i, module_id = rownames(temp),
        module_name = temp$`temp[[1]]$MODULE`
      ))
      if (ncol(temp) > 1) {
        return(temp)
      }
    }
    return(NULL)
  }))

  metabolite_modules <- metabolite_modules %>% left_join(temp2,
    by = "KEGG",
    relationship = "many-to-many"
  )

  # use lapply() to apply KEGGREST::keggGet() to each unique module ID in parallel
  temp2 <- do.call(rbind, lapply(unique(na.omit(metabolite_modules$module_id)), function(i) {
    temp <- tryCatch(
      {
        Sys.sleep(0.2) # add a delay of 1 second between requests
        KEGGREST::keggGet(i)
      },
      error = function(e) NULL
    )
    temp <- strsplit(temp[[1]][["CLASS"]], ";")
    temp <- as.data.frame(cbind(
      module_id = i, upper_hierarchy = str_squish(temp[[1]][1]),
      middle_hierarchy = str_squish(temp[[1]][2]),
      lower_hierarchy = str_squish(temp[[1]][3])
    ))
    return(temp)
  }))

  # join with previous data frame
  metabolite_modules <- metabolite_modules %>% left_join(temp2,
    by = "module_id",
    relationship = "many-to-many"
  )

  metabolite_modules <- metabolite_modules[!is.na(metabolite_modules$module_id), ]
  ORA_dataframes[["annotation"]] <- metabolite_modules

  if (update_background == TRUE) {
    message("You chose to retrieve the background KEGG information.
            This request takes some time, please be patient.")
    # background information of all metabolites in KEGG modules ####
    # annotate all metabolites to a module
    # get list of all KEGG modules
    all_modules <- as.data.frame(KEGGREST::keggList("module"))
    all_modules <- as.data.frame(cbind(
      module_id = rownames(all_modules),
      module_name = all_modules$`KEGGREST::keggList("module")`
    ))

    modules_compounds <- as.data.frame(matrix(ncol = 6))
    colnames(modules_compounds) <- c(
      "module_id", "module_name", "kegg_id",
      "upper_hierarchy", "middle_hierarchy", "lower_hierarchy"
    )

    # use lapply() to apply KEGGREST::keggGet() to each module ID in parallel
    modules_compounds <- do.call(rbind, lapply(all_modules$module_id, function(i) {
      temp <- tryCatch(
        {
          Sys.sleep(0.2) # add a delay of 1 second between requests
          KEGGREST::keggGet(i)
        },
        error = function(e) NULL
      )
      if (rownames(as.data.frame(temp[[1]][["ENTRY"]])) == "Pathway" &
        length(intersect(names(temp[[1]]), "COMPOUND")) > 0) {
        temp2 <- as.data.frame(temp[[1]][["COMPOUND"]])
        temp3 <- strsplit(temp[[1]][["CLASS"]], ";")
        temp3 <- as.data.frame(cbind(
          module_id = i, KEGG = rownames(temp2),
          module_name = temp[[1]][["NAME"]],
          upper_hierarchy = str_squish(temp3[[1]][1]),
          middle_hierarchy = str_squish(temp3[[1]][2]),
          lower_hierarchy = str_squish(temp3[[1]][3])
        ))
        return(temp3)
      }
      return(NULL)
    }))

    modules_compounds <- modules_compounds[!is.na(modules_compounds$module_id), ]
    ORA_dataframes[["background"]] <- modules_compounds
  }
  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["KEGG_annotations"]] <- ORA_dataframes
    return(data)
  } else {
    # otherwise, return the list of fits
    return(ORA_dataframes)
  }
}

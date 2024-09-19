#' Retrieve background and annotation information for over-representation analysis
#' (ORA)
#'
#' Uses the package KEGGREST to retrieve background and annotation information
#' needed for over-representation analysis. As KEGGREST only allows 10 queries
#' per second this might take some time to run, depending on the size of the
#' dataset and organism. The user should check afterwards if all functional modules are
#' applicable for the analysis of the dataset (p.e. organism, tissue).
#'
#' @param data data set to be analyzed with ORA. Must at least contain a column
#'                      with KEGG IDs
#' @param kegg column name of "data" that holds KEGG IDs
#' @param metabolite_name column name of "data" that holds metabolite names
#' @param update_background logical. Should the background information be updated?
#' If TRUE this may take some time.
#'
#' @seealso [ORA_hypergeometric()]
#'
#' @return a list with dataframes "background" and "annotation" needed for ORA
#' @export
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGREST keggList
#' @importFrom SummarizedExperiment colData
#' @importFrom stats na.omit
#' @import dplyr
#' @import stringr
#'
#' @examples
#' data("data_sim")
#' data_sim <- as.data.frame(SummarizedExperiment::colData(data_sim))
#' data <- data_sim[data_sim$metabolite == "ATP", ]
#' ORA_dataframes <- get_ORA_dataframes(
#'   data = data,
#'   kegg = "KEGG",
#'   metabolite_name = "metabolite",
#'   update_background = FALSE
#' )
#' head(ORA_dataframes[["annotation"]])
get_ORA_dataframes <- function(data, kegg = "KEGG",
                               metabolite_name = "metabolite",
                               update_background = FALSE) {
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    data <- as.data.frame(SummarizedExperiment::colData(data))
  }

  # bind variables to function
  ORA_dataframes <- list()
  metabolite_modules <- NULL
  temp2 <- NULL
  temp <- NULL
  temp3 <- NULL
  modules_compounds <- NULL


  # hand KEGG IDs to KEGGREST, KEGGREST can only do 10 queries at a time
  metabolite_modules <- unique(as.data.frame(cbind(
    metabolite = data[[metabolite_name]],
    KEGG = data[[kegg]]
  )))

  temp2 <- as.data.frame(cbind(KEGG = NA, module_id = NA, module_name = NA))
  for (i in unique(data[[kegg]])) {
    # only if there is a KEGG ID
    if (!is.na(i)) {
      cat(paste0(i," "))
      # if KEGGREST runs into warnings the loop will brake
      temp <- tryCatch(KEGGREST::keggGet(i), error=function(e) NULL)
      # only modules containing only human pathways
      temp <- as.data.frame(temp[[1]]$MODULE)
      # extract module id and name
      temp <- as.data.frame(cbind(
        KEGG = i, module_id = rownames(temp),
        module_name = temp$`temp[[1]]$MODULE`
      ))
      # only if there is a module entry to KEGG ID
      if (ncol(temp) > 1) {
        temp2 <- rbind(temp, temp2)
      }
    }
  }

  metabolite_modules <- metabolite_modules %>% left_join(temp2,
    by = "KEGG",
    relationship = "many-to-many"
  )
  rm(i, temp, temp2)

  # we also want the hierarchy of modules, so we hand the identified module ids
  # again to keggGet and extract the entries of $class

  temp2 <- as.data.frame(cbind(
    module_id = NA, upper_hierarchy = NA, middle_hierarchy = NA,
    lower_hierarchy = NA
  ))

  for (i in unique(na.omit(metabolite_modules$module_id))) {
    temp <- tryCatch(KEGGREST::keggGet(i), error=function(e) NULL)
    # only select modules that only contain human pathways
    # entry Class is semicolon seperated in order of hierarchy
    temp <- strsplit(temp[[1]][["CLASS"]], ";")
    # assign class entries to hierarchies
    temp <- as.data.frame(cbind(
      module_id = i, upper_hierarchy = str_squish(temp[[1]][1]),
      middle_hierarchy = str_squish(temp[[1]][2]),
      lower_hierarchy = str_squish(temp[[1]][3])
    ))
    # collect in one dataframe
    temp2 <- rbind(temp2, temp)
  }

  # join with previous dataframe
  metabolite_modules <- metabolite_modules %>% left_join(temp2,
    by = "module_id",
    relationship = "many-to-many"
  )
  rm(i, temp, temp2)

  metabolite_modules <- metabolite_modules[!is.na(metabolite_modules$module_id), ]
  ORA_dataframes[["annotation"]] <- metabolite_modules
  rm(metabolite_modules)

  if (update_background == TRUE) {
    # background information of all metabolites in KEGG modules ####
    # annotate all metabolites to a module
    # get list of all KEGG modules
    all_modules <- as.data.frame(KEGGREST::keggList("module"))
    all_modules <- as.data.frame(cbind(
      module_id = rownames(all_modules),
      module_name = all_modules$`keggList("module")`
    ))

    modules_compounds <- as.data.frame(matrix(ncol = 5))
    colnames(modules_compounds) <- c(
      "module_id", "kegg_id",
      "upper_hierarchy", "middle_hierarchy", "lower_hierarchy"
    )
    # KEGGREST only allows for 10 simultaneous queries
    for (i in all_modules$module_id)
    {
      # get kegg Entry of module
      temp <- tryCatch(KEGGREST::keggGet(i), error=function(e) NULL)
      # some module_ids contain information about complexes but not pathway modules
      # and only select modules where we get compound information
      if (rownames(as.data.frame(temp[[1]][["ENTRY"]])) == "Pathway" &
        length(intersect(names(temp[[1]]), "COMPOUND")) > 0) {
        # retrieve list of compounds (=metabolite) in module
        temp2 <- as.data.frame(temp[[1]][["COMPOUND"]])
        # retrieve hierachical information
        temp3 <- strsplit(temp[[1]][["CLASS"]], ";")
        # assign class entries to hierachies
        temp3 <- as.data.frame(cbind(
          module_id = i, kegg_id = rownames(temp2),
          upper_hierarchy = str_squish(temp3[[1]][1]),
          middle_hierarchy = str_squish(temp3[[1]][2]),
          lower_hierarchy = str_squish(temp3[[1]][3])
        ))

        # combine information on all modules in one dataframe
        modules_compounds <- rbind(modules_compounds, temp3)
      }
    }
    rm(i, temp, temp2, temp3)
    rm(all_modules)
    modules_compounds <- modules_compounds[!is.na(modules_compounds$module_id), ]
    ORA_dataframes[["background"]] <- modules_compounds
    rm(modules_compounds)
  }
  return(ORA_dataframes)
}

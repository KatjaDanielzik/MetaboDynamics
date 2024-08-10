# # We want to map metabolites to a coarse-grain version of pathways. In the KEGG
# # database we can accesss modules: accumulations of pathways with biological
# # function
#
# # experimental metabolites in modules ####
#
# # load package KEGGREST which allows an access to the KEGG API
# if (!require("KEGGREST", quietly = TRUE)) {
#   BiocManager::install("KEGGREST")
# }
# library("KEGGREST")
# library(stringr)
#
# # load list of metabolites from script (script/CAS_ID_mapping_for_graphite.R)
# # which is dependent on script/input_metaboanalyst.R and contains metabaolite
# # names and KEGG IDS
# name_map_HMDB_CAS <- readr::read_csv("/media/home/5-1-package pipeline/data/name_map_HMDB_CAS.csv")
#
# # hand KEGG IDs to KEGGREST, KEGGREST can only do 10 queries at a time
# metabolite_modules <- as.data.frame(cbind(
#   metabolite = name_map_HMDB_CAS$Query,
#   KEGG = name_map_HMDB_CAS$KEGG
# ))
# # get list of human pathways to only select modules that contain human pathways
# human_pathways <- as.data.frame(keggList("pathway", "hsa"))
# human_pathways <- as.data.frame(cbind(
#   pathway_id = rownames(human_pathways),
#   pathway_name = human_pathways$`keggList("pathway", "hsa")`
# ))
# ## pathway ids are represented as mapID in the list of modules -> change "hsa" to "map"
# human_pathways$pathway_id <- gsub("hsa", "map", human_pathways$pathway_id)
#
# temp2 <- as.data.frame(cbind(KEGG = NA, module_id = NA, module_name = NA))
# for (i in name_map_HMDB_CAS$KEGG) {
#   # only if there is a KEGG ID
#   if (!is.na(i)) {
#     temp <- keggGet(i)
#     # only modules containing only human pathways
#     temp <- as.data.frame(temp[[1]]$MODULE)
#     # extract module id and name
#     temp <- as.data.frame(cbind(
#       KEGG = i, module_id = rownames(temp),
#       module_name = temp$`temp[[1]]$MODULE`
#     ))
#     # only if there is a module entry to KEGG ID
#     if (ncol(temp) > 1) {
#       temp2 <- rbind(temp, temp2)
#     }
#   }
# }
#
# library(dplyr)
# metabolite_modules <- metabolite_modules %>% left_join(temp2, by = "KEGG")
# rm(i, temp, temp2)
#
# # we also want the hierachy of modules, so we hand the identified module ids
# # again to keggGet and extract the entries of $class
#
# temp2 <- as.data.frame(cbind(
#   module_id = NA, upper_hierachy = NA, middle_hierachy = NA,
#   lower_hierachy = NA
# ))
#
#
# for (i in unique(na.omit(metabolite_modules$module_id))) {
#   temp <- keggGet(i)
#   # only select modules that only contain human pathways
#   # entry Class is semicolon seperated in order of hierachy
#   temp <- strsplit(temp[[1]][["CLASS"]], ";")
#   # assign class entries to hierachies
#   temp <- as.data.frame(cbind(
#     module_id = i, upper_hierachy = str_squish(temp[[1]][1]),
#     middle_hierachy = str_squish(temp[[1]][2]),
#     lower_hierachy = str_squish(temp[[1]][3])
#   ))
#   # collect in one dataframe
#   temp2 <- rbind(temp2, temp)
# }
#
# # join with previous dataframe
# metabolite_modules <- metabolite_modules %>% left_join(temp2, by = "module_id")
# rm(i, temp, temp2)
#
# # we are working in a human lung cancer cell line
# # "Xenobiotics biodegradation" physiological only present in liver
# metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Xenobiotics biodegradation"), ]
# metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Biosynthesis of other secondary metabolites"), ]
# metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Biosynthesis of terpenoids and polyketides"), ]

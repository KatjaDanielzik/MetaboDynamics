\name{metabolite_modules}
\alias{metabolite_modules}
\docType{data}
\title{
KEGG Query Results of experimental metabolites 
}
\description{
Using the package KEGGREST (https://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html) all experimental metabolites (see data("intra")) were queried with there KEGG-IDs and all functional modules recorded to which the metabolite is annotated in the KEGG-database.
Restricted to HUMAN pathways and DOES
NOT include functional modules "Xenobiotics biodegradation", "Biosynthesis of other secondary metabolites",
"Biosynthesis of terpenoids and polyketides"
}
\usage{data("metabolite_modules")}
\format{
  A data frame with 348 observations on the following 8 variables.
  \describe{
    \item{\code{...1}}{row number of the dataframe}
    \item{\code{metabolite}}{name of the experimental metabolite}
    \item{\code{KEGG}}{KEGG ID of the experimental metabolite}
    \item{\code{module_id}}{ID of the KEGG module to which the metabolite is annotated}
    \item{\code{module_name}}{name of the KEGG module to which the metabolite is annotated}
    \item{\code{upper_hierarchy}}{name of the highest hierachy level of module organisation}
    \item{\code{middle_hierarchy}}{name of the middle hierachy = functional module, p.e. "Amino acid metabolism"}
    \item{\code{lower_hierarchy}}{name of the lowest level of modules, this often contain only a couple pathways p.e. "Arginine and proline metabolism" }
  }
}
\source{
https://www.genome.jp/kegg/module.html
data <- get_ORA_annotations(longitudinalMetabolomics,update_background=TRUE)
metabolite_modules <- get_ORA_annotations["annotations"]

# we are working in a human lung cancer cell line
# "Xenobiotics biodegradation" physiological only present in liver
metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Xenobiotics biodegradation"), ]
metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Biosynthesis of other secondary metabolites"), ]
metabolite_modules <- metabolite_modules[-which(metabolite_modules$middle_hierachy == "Biosynthesis of terpenoids and polyketides"), ]

}
\seealso{
\link{modules_compounds}/ 
\code{\link[=get_ORA_annotations]{get_ORA_annotations()}}
}

\keyword{datasets}

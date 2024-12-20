\name{modules_compounds}
\alias{modules_compounds}
\docType{data}
\title{
Background KEGG Query Results Of Functional Modules
}
\description{
Using the package KEGGREST (https://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html) a list of all
KEGG-modules (KeggList("module")) including their upper, middle and lower hierachy as given by the KEGG-database and the corresponding annotated metabolites was queried.
}
\usage{data("modules_compounds")}
\format{
  A data frame with 1988 observations on the following 6 variables.
  \describe{
    \item{\code{KEGG}}{KEGG ID of a metabolite annotated to a functional module}
    \item{\code{upper_hierarchy}}{name of the highest hierachy level of module organisation}
    \item{\code{middle_hierarchy}}{name of the middle hierachy = functional module, p.e. "Amino acid metabolism"}
    \item{\code{lower_hierarchy}}{name of the lowest level of modules, this often contain only a couple pathways p.e. "Arginine and proline metabolism" }
    \item{\code{module_id}}{the ID of the KEGG functional module}
    \item{\code{module_name}}{name of the KEGG module}
  }
}
\source{
https://www.genome.jp/kegg/module.html
data <- get_ORA_annotations(longitudinalMetabolomics,update_background=TRUE)
modules_compunds <- get_ORA_annotations["background"]
}
\seealso{
\link{metabolite_modules}/ 
\code{\link[=get_ORA_annotations]{get_ORA_annotations()}}
}
\keyword{datasets}

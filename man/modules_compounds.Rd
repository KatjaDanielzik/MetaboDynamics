\name{modules_compounds}
\alias{modules_compounds}
\docType{data}
\title{
Background KEGG Query Results Of Functional Modules
}
\description{
Using the package KEGGREST (https://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html) a list of all
KEGG-modules (KeggList("module")) including their upper, middle and lower hierachy as given by the KEGG-database and the corresponding annotated metabolites was queried.
Restricted to HUMAN pathways and DOES
NOT include functional modules "Xenobiotics biodegradation", "Biosynthesis of other secondary metabolites",
"Biosynthesis of terpenoids and polyketides"
}
\usage{data("modules_compounds")}
\format{
  A data frame with 1988 observations on the following 6 variables.
  \describe{
    \item{\code{...1}}{is a running ID of the dataframe row}
    \item{\code{module_id}}{the ID of the KEGG functional module}
    \item{\code{kegg_id}}{KEGG ID of a metabolite annotated to a functional module}
    \item{\code{upper_hierarchy}}{name of the highest hierachy level of module organisation}
    \item{\code{middle_hierarchy}}{name of the middle hierachy = functional module, p.e. "Amino acid metabolism"}
    \item{\code{lower_hierarchy}}{name of the lowest level of modules, this often contain only a couple pathways p.e. "Arginine and proline metabolism" }
  }
}
\details{
script to obtain data is in R/modules_compounds.R (last executed September 2023)
}
\source{
https://www.genome.jp/kegg/module.html
}
\keyword{datasets}

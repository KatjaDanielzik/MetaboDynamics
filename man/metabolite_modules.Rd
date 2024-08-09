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
    \item{\code{upper_hierachy}}{name of the highest hierachy level of module organisation}
    \item{\code{middle_hierachy}}{name of the middle hierachy = functional module, p.e. "Amino acid metabolism"}
    \item{\code{lower_hierachy}}{name of the lowest level of modules, this often contain only a couple pathways p.e. "Arginine and proline metabolism" }
  }
}
\details{
script to obtain data is in R/metabolites_modules.R (last executed September 2023)
}
\source{
https://www.genome.jp/kegg/module.html
}
\examples{
data("metabolite_modules")
unique(metabolite_modules[metabolite_modules$metabolite=="ATP",]$module_name)
}
\keyword{datasets}

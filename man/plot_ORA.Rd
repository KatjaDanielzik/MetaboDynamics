% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ORA.R
\name{plot_ORA}
\alias{plot_ORA}
\title{Plot results of over-representation analysis with ORA_hypergeometric()}
\usage{
plot_ORA(data, tested_column = "middle_hierarchy")
}
\arguments{
\item{data}{result dataframe from ORA_hypergeometric() or \link[SummarizedExperiment]{SummarizedExperiment}
object where the ORA_hypergeometric() results are stored in metadata(data)
under "ORA_tested_column"}

\item{tested_column}{KEGG module hierarchy level on which ORA was executed}
}
\value{
a plot of the over-representation analysis
}
\description{
Plot results of over-representation analysis with ORA_hypergeometric()
}
\examples{
data("longitudinalMetabolomics")
data("modules_compounds")
head(modules_compounds)
data("metabolite_modules")
head(metabolite_modules)
# middly hierachy
longitudinalMetabolomics <- ORA_hypergeometric(
  data = longitudinalMetabolomics,
  annotations = metabolite_modules,
  background = modules_compounds,
  tested_column = "middle_hierarchy"
)
plot_ORA(longitudinalMetabolomics)
}
\seealso{
do over-represenation analysis of KEGG functional modules \code{\link[=ORA_hypergeometric]{ORA_hypergeometric()}}
}

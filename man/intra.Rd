\name{intra}
\alias{intra}
\docType{data}
\title{
Intracellular Metabolite Concentrations
}
\description{
This dataset contains intracellular metabolite concentrations normalized by cell counts in three replicates (each in one cell well) each
of three radiation doses (0,2,10Gy) at four timepoints (0h=before iradiation, 1h, 6h, 24h post-irradiation).
Untargeted LC-MS was performed to obtain metabolite concentrations for every sample.
The experiment was conducted by Dr. Johann Matschke and LC-MS was performed by Prof. Dr. Thomer Shlomi. Maven was used as a software for peak annotation. We consider peakintensity/cell count = metabolite concentration per cell. Samples were handled in paarallel by the same experimentor.
}
\usage{data("intra")}
\format{
  A data frame with 2940 observations on the following 7 variables.
  \describe{
    \item{\code{dose}}{a character vector specifying the radiation dose}
    \item{\code{time}}{a character vector specifying the time post-irradiation. 0h refers to directly before irradiation}
    \item{\code{sample.ID}}{specifies the cell well used for the experiment (one per replicate,timepoint and dose)}
    \item{\code{cells}}{number of cells in the cell well}
    \item{\code{metabolite}}{Name of the metabolite as annotated by MAVEN}
    \item{\code{value}}{measured peak intensity during LC-MS}
    \item{\code{cpc}}{Concentration Per Cell= Value/cells}
  }
}
%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{
Johann Matschke
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data("intra")
head(intra)
}
\keyword{datasets}

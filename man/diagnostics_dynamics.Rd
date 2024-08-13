\name{diagnostics_dynamics}
\alias{diagnostics_dynamics}
\docType{data}
\title{
Diagnostics of fits that come from function fits_dynamics_model().
}
\description{
Holds the diagnostic criteria: 1) number of divergent transitions("divergences"),
exceeding maximum treedepth ("treedepth_error"), rhat and number of effective
samples (neff).
}
\usage{data("diagnostics_dynamics")}
\format{
  The format is:
List of 1
 $ model_diagnostics:'data.frame':	294 obs. of  12 variables:
  ..$ metabolite.ID
  ..$ condition
  ..$ divergences
  ..$ treedepth_error
  ..$ rhat_mu1_mean
  ..$ rhat_mu2_mean
  ..$ rhat_mu3_mean
  ..$ rhat_mu4_mean
  ..$ neff_mu1_mean
  ..$ neff_mu2_mean
  ..$ neff_mu3_mean
  ..$ neff_mu4_mean
}
\source{
/R/prepare_dynamic_model.R
}
\examples{
data(diagnostics_dynamics)
## maybe str(diagnostics_dynamics) ; plot(diagnostics_dynamics) ...
}
\keyword{datasets}

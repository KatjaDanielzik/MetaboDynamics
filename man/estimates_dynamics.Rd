\name{estimates_dynamics}
\alias{estimates_dynamics}
\docType{data}
\title{
Diagnostics of fits that come from function fits_dynamics_model().
}
\description{
Holds for every experimental condition a dataframe collecting 1) mu_mean at
every time point. 2) Highest density interval ("...lower","...upper") of mu,
sigma and lambda (see /inst/stan/m_ANOVA_partial_pooling.stan) as well as 3)
100 samples from the posterior of mu_means at every timepoint and 4) differences
between two consecutive timepoints mu2_mean-mu1_mean ("delta...")
}
\usage{data("estimates_dynamics")}
\format{
  The format is:
List of 3
 $ 0Gy :'data.frame':	98 obs. of  439 variables:
  ..$ condition
  ..$ metabolite.ID
  ..$ metabolite
  ..$ mu1_mean
  ..$ mu2_mean
  ..$ mu3_mean
  ..$ mu4_mean
  ..$ mu1_lower
  ..$ mu2_lower
  ..$ mu3_lower
  ..$ mu4_lower
  ..$ mu1_higher
  ..$ mu2_higher
  ..$ mu3_higher
  ..$ mu4_higher
  ..$ sigma1_mean
  ..$ sigma2_mean
  ..$ sigma3_mean
  ..$ sigma4_mean
  ..$ sigma1_lower
  ..$ sigma2_lower
  ..$ sigma3_lower
  ..$ sigma4_lower
  ..$ sigma1_higher
  ..$ sigma2_higher
  ..$ sigma3_higher
  ..$ sigma4_higher
  ..$ lambda_mean
  ..$ lambda_lower
  ..$ lambda_higher
  ..$ delta12_mean
  ..$ delta12_lower
  ..$ delta12_higher
  ..$ delta23_mean
  ..$ delta23_lower
  ..$ delta23_higher
  ..$ delta34_mean
  ..$ delta34_lower
  ..$ delta34_higher
  ..$ mu1_sample1
  ..$ mu2_sample1
  ..$ mu3_sample1
  ..$ mu4_sample1
 $ 2Gy :
 $ 10Gy:
}
\source{
/R/prepare_dynamic_model.R
}
\examples{
data(estimates_dynamics)
## maybe str(estimates_dynamics) ; plot(estimates_dynamics) ...
}
\keyword{datasets}

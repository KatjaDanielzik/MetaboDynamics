#' extract_estimates_dynamics
#'
#' Extracts the mean concentrations (mu) at every timepoint from the dynamic model fit as well as the 95% HDI. As well as the estimated standard deviation of metabolite concentrations at every time point (sigma) and the pooled standard deviation of every metabolite over all timepoints (lambda).
#' Additionally samples from the posterior of mu can be drawn. This can be helpful if p.e. one wants to estimate the clustering precision. Lambda can be used for clustering algorithms such as VSClust that also take the variance into account.
#'
#' @param data dataframe used for modeling
#' @param M number of metabolites, default requires a column in data named "metabolite"
#' @param t number of unique timepoints in data, the default require a column named "time"
#' @param condition name of column in dataframe data that specifies the experimental condition, default requires a column named "dose"
#' @param fits lits of model fits for which estimates should be extracted
#' @param iter how many iterations were used to fit the dynamic model
#' @param warmup how many warm-up iterations were used to fit the dynamic model
#' @param chains how many chains were used to fit the dynamic model
#' @param samples how many posterior samples should be drawn (p.e. for check of clustering precision)
#'
#' @import dplyr
#' @importFrom stats runif
#'
#' @return a list of dataframes (one per experimental condition) that contains the estimates at the timepoints and samples from the posterior (number as specified in samples)
#' @export
#'
#' @examples
#' #' data("intra")
#' # only run after fit_dynamics_model(data = intra[intra$dose=="0Gy",],
## # cpc = "log_cpc_stand", condition = "dose", max_treedepth = 14, adapt_delta = 0.999, iter = 4000, cores = 7): see Vignette and documentation
#' #of function
#' # extract_estimates(data=intra,fits=fits_dynamics,iter=4000)
extract_estimates_dynamics<-function(data,M=length(unique(data$metabolite)),t=length(unique(data$time)),condition="dose",fits,iter=2000,warmup=iter/4,chains=4,samples=100){

  conditions <- unique(data[[condition]])

  dynamics <- list()

  # get estimates
  for (i in names(fits))
  {
    # generate dataframe for storage
    dynamics_log_cpc <- as.data.frame(cbind(condition=i,metabolite.ID=1:M))
    # select only models for which delta mu with 0Gy/0h was computed
    fit <- fits[[i]]
    # generate n=samples random draw numbers from posterior, but same draw
    # from every delta alpha as they are dependent on each other
    x <- floor(runif(samples, min=1, max=(iter-warmup)*chains))

    # posterior summary
    pS <- as.data.frame(rstan::summary(fit)$summary)

    # extract posterior samples
    mu_posterior <- as.matrix(fit)


    for (m in 1:M){
      # single loops to get reasonable order for clustering
      for (j in 1:t){
        # get means of estimated means at timepoints
        dynamics_log_cpc[m,paste0("mu",j,".mean")] <- pS[paste0("mu[",m,",",j,"]"),"mean"]
      }
      for(j in 1:t){
        # lower border of 95% CrI
        dynamics_log_cpc[m,paste0("mu",j,"_lower")] <- pS[paste0("mu[",m,",",j,"]"),"2.5%"]
      }
      for (j in 1:t){
        # higher border of 95% CrI
        dynamics_log_cpc[m,paste0("mu",j,"_higher")] <- pS[paste0("mu[",m,",",j,"]"),"97.5%"]
      }
      for (j in 1:t){
        # get means of sigma
        dynamics_log_cpc[m,paste0("sigma",j,"_mean")] <- pS[paste0("sigma[",m,",",j,"]"),"mean"]
      }
      for(j in 1:t){
        # lower border of 95% CrI
        dynamics_log_cpc[m,paste0("sigma",j,"_lower")] <- pS[paste0("sigma[",m,",",j,"]"),"2.5%"]
      }
      for (j in 1:t){
        # higher border of 95% CrI
        dynamics_log_cpc[m,paste0("sigma",j,"_higher")] <- pS[paste0("sigma[",m,",",j,"]"),"97.5%"]
      }
      # get means of lambda
      dynamics_log_cpc[m,"lambda_mean"] <- pS[paste0("lambda[",m,"]"),"mean"]
      # lower border of 95% CrI
      dynamics_log_cpc[m,"lambda_lower"] <- pS[paste0("lambda[",m,"]"),"2.5%"]
      # higher border of 95% CrI
      dynamics_log_cpc[m,"lambda_higher"] <- pS[paste0("lambda[",m,"]"),"97.5%"]

      # delta mus
      # we have one less delta_mu than timepoints
      for (j in 1:(t-1))
      {
        # mean
        dynamics_log_cpc[m,paste0("delta",j,j+1,"_mean")] <- pS[paste0("delta_mu[",m,",",j,"]"),"mean"]
        # lower border of 95% CrI
        dynamics_log_cpc[m,paste0("delta",j,j+1,"_lower")] <- pS[paste0("delta_mu[",m,",",j,"]"),"2.5%"]
        # higher border of 95% CrI
        dynamics_log_cpc[m,paste0("delta",j,j+1,"_higher")] <- pS[paste0("delta_mu[",m,",",j,"]"),"97.5%"]
      }


      for (k in 1:length(x))
      {
        for(j in 1:t){
          dynamics_log_cpc[m,paste0("mu",j,"_sample",k)] <- mu_posterior[,paste0("mu[",m,",",j,"]")][x[k]]
        }
      }
    }

    # get correct assignment of metabolites and metabolite.ID used in modelling
    metabolites <- as.data.frame(cbind(metabolite=unique(intra$metabolite),metabolite.ID=as.numeric(as.factor(intra$metabolite))))

    dynamics_log_cpc <- left_join(dynamics_log_cpc,metabolites,by="metabolite.ID")
    dynamics_log_cpc <- dynamics_log_cpc%>%relocate("metabolite",.after="metabolite.ID")
    dynamics_log_cpc <- unique(dynamics_log_cpc)
    dynamics[[i]] <- dynamics_log_cpc
  }
  return(dynamics)
  rm(pS,fit,i,mu_posterior,j,k,x,conditions,metabolites,dynamics_loc_cpc)
}

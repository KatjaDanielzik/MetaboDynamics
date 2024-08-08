#' Extracts diagnostic criteria from dynamics model
#'
#' gathers number of divergences, rhat values, number of effectives samples (n_eff) and provides plots for diagnostics criteria as well ad Posterior predictive checks
#'
#' @param data dataframe used to fit dynamics model
#' @param N number of rows in used dataframe
#' @param M number of metabolites in experimental dataset
#' @param t number of timepoints in experimental dataset
#' @param max_treedepth max_treedepth setting used to fit dynamic model
#' @param adapt_delta adapt_delta setting used to fit dynamic model
#' @param fits list of models for which diagnostics should be extracted, is the object that gets returned by fit_dynamics_model()
#' @param warmup number of warmup iterations used for model fit
#' @param iter number of iterations used for model fit
#' @param chains number of chains used for model fit
#'
#' @return a list which contains diagnostics criteria of all conditions in a dataframe (named "model_diagnostics") and one dataframe per condition that contains necessary information for Posterior predictive check (named "PPC_condition"). Additionally plots for diagnostics and PPC named "plot_criteria" and "plot_PPC_condition"
#' @export
#'
#' @examples
#' extract_diagnostics_dynamics(data=intra)
#'

extract_diagnostics_dynamics <- function(data,N=nrow(data),M=length(unique(data$metabolite)),t=length(unique(data$time)),max_treedepth=10,adapt_delta=0.95,iter=2000,warmup=iter/4,chains=4,fits="fits_dynamcis"){
  loadNamespace("ggplot2")
  loadNamespace("tidyr")
  # create list to store all subsequent results
  list_diagnostics <- list()

  # diagnostic criteria from model
  diagnostics_dynamics <- as.data.frame(matrix(ncol=13))
  colnames(diagnostics_dynamics) <- c("label","treedepth","adapt_delta","divergences","treedepth_error","rhat_mu1_mean","rhat_mu2_mean","rhat_mu3_mean","rhat_mu4_mean","n_eff_mu1_mean","n_eff_mu2_mean","n_eff_mu3_mean","n_eff_mu4_mean")

  for (i in names(fits)){
    # temporary file
    fit <- fits[[i]]
    # gather diagnostic criteria
    divergences <- rstan::get_num_divergent(fit)
    max_treedepth <- rstan::get_num_max_treedepth(fit)
    rhat <- rstan::summary(fit)$summary[,"Rhat"]
    n_eff <- rstan::summary(fit)$summary[,"n_eff"]
    label <- i

    # gather in one dataframe
    for (m in 1:M){
      start <- (m-1)*t+1
      end <- (m*t)
      temp <- cbind(label=label,treedepth=max_treedepth,adapt_delta=adapt_delta,divergences=divergences, treedepth_error=max_treedepth,rhat_mu1_mean=rhat[start],rhat_mu2_mean=rhat[start+1],rhat_mu3_mean=rhat[start+2],rhat_mu4_mean=rhat[end],n_eff_mu1_mean=n_eff[start],n_eff_mu2_mean=n_eff[start+1],n_eff_mu3_mean=n_eff[start+2],n_eff_mu4_mean=n_eff[end])
      diagnostics_dynamics <- rbind(diagnostics_dynamics,temp)
    }

    # Posterior predictive check
    # turn y_rep from model fit into long format and add metabolite.ID and time.ID
    posterior <- as.data.frame(fit,pars="y_rep")
    posterior <- pivot_longer(posterior,names_to="parameter",values_to = "posterior",cols = 1:all_of(M*t))
    posterior$metabolite.ID <- as.numeric(rep(rep(1:M,t),(iter-warmup)*chains))
    posterior$time.ID <- as.factor(rep(rep(1:t,each=M),(iter-warmup)*chains))

    # assign metabolite and time id to data
    PPC <- data
    PPC$metabolite.ID <- as.numeric(as.factor(PPC$metabolite)) # equivalent to model coding
    PPC$time.ID <- as.factor(as.numeric(as.factor(as.numeric(PPC$time))))

    list_diagnostics[[paste0("PCC_",i)]]<-PPC
    # visualize PPC

  }

  list_diagnostics[["model_diagnostics"]] <- diagnostics_dynamics
  #cleanup
  rm(i,fit,m,label,n_eff,max_treedepth,divergences,rhat,PPC)

  # visualize
  list_diagnostics[["plot_divergences"]]<-
    ggplot2::ggplot(diagnostics_dynamics, aes(x=label,y=as.numeric(divergences)))+
    geom_violin()+
    xlab("condition")+
    theme_bw()+
    ylab("number of divergent transitions")+
    ggtitle("diagnostics_dynamics","divergent transitions?")
  list_diagnostics[["plot_treedepth_error"]]<-
  ggplot2::ggplot(diagnostics_dynamics, aes(x=label,y=as.numeric(treedepth_error)))+
    geom_violin()+
    xlab("condition")+
    theme_bw()+
    ylab("number of exceeded treedepth")+
    ggtitle("diagnostics_dynamics","maximum treedepth exceeded?")


  temp <- diagnostics_dynamics %>% pivot_longer(cols=6:(6-1+t),names_to = "rhat.mu" ,values_to ="rhat.value")
  list_diagnostics[["plot_rhat"]]<-
  ggplot2::ggplot(temp, aes(x=as.numeric(as.factor(rhat.mu)),y=as.numeric(rhat.value)))+
    geom_violin()+
    theme_bw()+
    facet_wrap(~label)+
    ylab("Rhat")+
    xlab("timepoint")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_hline(yintercept=1.01,col="red")+
    ggtitle("diagnostics_dynamics","Rhat < 1.01 ?")

  # number of effective samples
  temp <- diagnostics_dynamics %>% pivot_longer(cols=10:(10-1+t),names_to = "neff.mu" ,values_to ="neff.value" )
  list_diagnostics[["plot_neff"]]<-
  ggplot2::ggplot(temp, aes(x=as.numeric(as.factor(neff.mu)),y=as.numeric(neff.value)))+
    geom_violin()+
    theme_bw()+
    facet_wrap(~label)+
    ylab("number of effective samples")+
    xlab("timepoint")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_hline(yintercept=100,col="red")+
    ggtitle("diagnostics_dynamics","number of effective samples >100 ?")

  return(list_diagnostics)

    # cleanup
  rm(temp)
}


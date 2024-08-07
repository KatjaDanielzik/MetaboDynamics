
#' Fits dynamic model
#'
#' MetaboDynamics provides a hierachical model for robust estimation of mean
#' concentrations over time of single metabolites at single experimental conditions.
#'
#' @param data concentration table containing the columns "metabolite", "condition", "cpc_stand" for default
#' @param metabolite column of "data" that contains the metabolite names or IDs
#' @param time column of "time" that contains time as numeric, make sure your time column is ordered from lowest to highest for the model to work
#' @param condition column of "data" that contains the experimental conditions
#' @param cpc column of "data" that contains the concentrations per cell, centered and normalized per metabolite and experimental condition (mean=0, sd=1)
#'
#' @return returns one modelfit named fit_condition per experimental condition
#' @export
#'
#' @examples
#' fit_dynamics_model(intra)
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MetaboDynamics

fit_dynamics_model <- function(data=intra,metabolite="metabolite",time="time",condition="dose",cpc="log_cpc_stand"){
# get unique experimental conditions
conditions <- unique(data[[condition]])
fits <- list()
# loop over conditions
for(i in 1:length(conditions)){
  # subset dataframe
temp <- data[data[[condition]]==conditions[i],]

# fit model
fit <- rstan::sampling(object=stanmodels$m_ANOVA_partial_pooling,
                       data=list(y=temp[[cpc]],
                                 t=length(unique(temp[[time]])),
                                 M=length(unique(temp[[metabolite]])),
                                 N=nrow(temp),
                                 # Vector of metabolite IDs
                                 Me=as.numeric(as.factor(temp[[metabolite]])),
                                 # Vector indicating which row belongs to which timestep
                                 X=as.numeric(as.factor(as.numeric(temp[[time]])))),
                       chains = 4,
                       iter = 4000,
                      # increase warmup so that algorithm probably chooses
                      #smaller step sizes for sampling the posterior
                       warmup = 1000,
                       algorithm = "NUTS",
                       cores = 4,
                       # increase adapt_delta (target average proposal acceptance probability)
                       # to decrease step size
                       control=list(adapt_delta=0.999,max_treedepth=14))
# assign condition name to fit and store in list
fits[[conditions[i]]] <- assign(paste0("fit","_",conditions[i]),fit)
}
# return list of fits
return(fits)
#cleanup
rm(i,temp,conditions,fit,fits)
}




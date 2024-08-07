#' Title
#'
#' @param data
#' @param metabolite
#' @param condition
#'
#' @return
#' @export
#'
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib MetaboDynamics
#'
#' @examples
fit_dynamics_model <- function(data,metabolite,condition){
rstan::sampling(stanmodels$m_ANOVA_partial_pooling)
}

# N <- nrow(temp)
# t <- length(unique(temp$time))
# M <- length(unique(temp$metabolite))
# r <- 3 # number of replicates
# fit_complete_pooling <- rstan::sampling(object = compiled_model_complete_pooling,
#                                         data = list(y = temp$log_cpc_stand,
#                                                     t = length(unique(temp$time)),
#                                                     M = length(unique(temp$metabolite)),
#                                                     N = nrow(temp),
#                                                     Me= as.numeric(as.factor(temp$metabolite)), # 98 metabolites * 4 time points
#                                                     # Vector indicating which row belongs to which timestep
#                                                     X=as.numeric(as.factor(as.numeric(temp$time)))),
#                                         chains = 4,
#                                         iter = 4000,
#                                         # increase warmup so that algorithm probably chooses
#                                         # smaller step sizes for sampling the posterior
#                                         warmup = 1000,
#                                         algorithm = "NUTS",
#                                         cores = 4,
#                                         # increase adapt_delta (target average proposal acceptance probability)
#                                         # to decrease step size
#                                         control=list(adapt_delta=.999,max_treedepth=14))

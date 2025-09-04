#' Extracts parameter estimates from numeric fit of Bayesian model of dynamics
#'
#' Extracts the mean concentrations (mu) at every time point from the dynamics model fit, the 95% highest density interval (HDI), the estimated standard deviation of metabolite concentrations at every time point (sigma), and the pooled standard deviation of every metabolite over all timepoints (lambda).
#' Additionally samples from the posterior of mu can be drawn. This can be helpful if p.e. one wants to estimate the clustering precision. Lambda can be used for clustering algorithms such as VSClust that also take the variance into account.
#'
#' @param data data frame or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used used to fit dynamics model, must contain a column specifying KEGG IDs, column named "condition" specifiyng the experimental condition and a column named "time" specifying the timepoints.
#' If it is a SummarizedExperiment object the dynamic fits must be stores in metadata(data)
#' under "dynamic_fits"
#' @param assay of the SummarizedExperiment object that was used to fit the dynamics
#' model
#' @param kegg column in "data" that contains the KEGG IDs or other identifier of metabolites
#' @param condition name of column in dataframe data that specifies the experimental condition
#' @param time column in "data" that contains the time point identifiers
#' @param metabolite column of "data" that contains the metabolite names or IDs
#' @param fit model fit for which estimates should be extracted
#'
#' @seealso Fit the dynamic model [fit_dynamics_model()].
#' Diagnostics of the dynamic model [diagnostics_dynamics()]
#' Visualization of estimates with [plot_estimates()]
#'
#' @import dplyr
#' @importFrom stats runif
#' @import SummarizedExperiment
#' @importFrom S4Vectors metadata
#'
#' @return a list of dataframes (one per parameters mu, sigma, lambda, delta_mu and euclidean distance) that contains
#' the estimates:
#' mu: is the estimated mean metabolite abundance
#' sigma: the estimated standard deviation of metabolite abundance
#' lambda: pooled sigma per condition
#' delta_mu: differences of mu between time points
#' euclidean_distance: estimated euclidean distance of time vectors of one metabolite between conditions
#' If data is a summarizedExperiment object the estimates are stored in
#' metadata(data) under "estimates_dynamics"
#' @export
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'                                    longitudinalMetabolomics$metabolite == "ATP"]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#'data <- estimates_dynamics(
#'   data = data
#' )
#' S4Vectors::metadata(data)[["estimates_dynamics"]]
estimates_dynamics <- function(data, assay = "scaled_log",
                               kegg = "KEGG", condition = "condition", time = "time",
                               metabolite = "metabolite",
                               fit = metadata(data)[["dynamic_fit"]]) {
  
  
  # Input checks
  if (!is.data.frame(data) & !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or colData of a SummarizedExperiment object")
  }
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    t <- nrow(rowData(data))
    data_df <- as.data.frame(cbind(
      as.data.frame(t(assays(data)[[assay]])),
      as.data.frame(colData(data))
    ))
    data_df <- data_df %>% pivot_longer(
      cols = seq_len(t), names_to = "time",
      values_to = "scaled_measurement"
    )
    fit <- metadata(data)[["dynamic_fit"]]
  }
  
  # convert potential tibbles into data frame
  if (is(data, "tbl")) {
    data <- as.data.frame(data)
  }
  if (is(data, "data.frame")) {
    data_df <- data
  }
  # check if all elements of fits are stanfit objects
  if (!inherits(fit, "stanfit")) {
    stop("'fit' must be a stanfit objects")
  }
  if (!all(vapply(list(time, condition, kegg), is.character, logical(1)))) {
    stop("'time', 'kegg' and 'condition' must be a character vector specifying a column name of data")
  }
  if (!all(c(time, kegg, condition) %in% colnames(data_df))) {
    stop("'data' must contain columns named 'condition', 'kegg' and 'time'")
  }

  # get number of metabolites, time points and conditions
  M <- length(unique(data_df[[metabolite]]))
  t <- length(unique(data_df[[time]]))
  C <- length(unique(data[[condition]]))


  estimates_data <- data.frame(
    metabolite = rep(rep(unique(data_df$metabolite),each=C),each=t),
    time = rep(rep(unique(data_df$time),each=C),M),
    condition = rep(rep(unique(data_df$condition),t),M)
  )
  
  # extract for mu
  mu <- rstan::summary(fit,pars="mu")$summary
  mu <- cbind(estimates_data, parameter = "mu",mu[,c("mean","2.5%","97.5%")])

  # extract for sigma
  sigma <- rstan::summary(fit,pars="sigma")$summary
  sigma <- cbind(estimates_data, parameter = "sigma", sigma[,c("mean","2.5%","97.5%")])
  
  # extract for lambda
  ## lambda only one per condition -> adapt estimates_data
  lambda_data <- data.frame(
    metabolite = rep(unique(data_df$metabolite),each=C),
    condition = rep(unique(data_df$condition),M)
  )

  lambda <- rstan::summary(fit,pars="lambda")$summary
  lambda <- cbind(lambda_data, parameter = "lambda", lambda[,c("mean","2.5%","97.5%")])
  
  # extract euclidean distances
  ## get possible dose combinations
  if(C>1&t>1){
  combinations <- t(combn(unique(data_df$condition), 2))

  distances_data <- data.frame(
    metabolite = rep(unique(data_df$metabolite),each=nrow(combinations)),
    condition_1 = combinations[,1],
    condition_2 = combinations[,2]
  )
  
  distances <- as.data.frame(rstan::summary(fit,pars="euclidean_distance")$summary)
  distances <- distances%>%na.omit() # posterior contains estiamtes that do not fullfull condition_1<conditon_2 due to stan technicalities
  distances <- cbind(distances_data,parameter="euclidean_distance",distances[,c("mean","2.5%","97.5%")])
  }
  if(t==1|C==1){
    distances <- "only possible for more than one time point and more than one condition"
  }
  
  # extract for delta_mu
  if (t > 1) {
    combinations <- t(combn(unique(data_df$time), 2))

    delta_mu_data <- data.frame(
      metabolite = rep(rep(unique(data_df$metabolite),each=nrow(combinations)),each=C),
      condition = rep(rep(unique(data_df$condition),each=nrow(combinations)),M),
      timepoint_1 = combinations[,1],
      timepoint_2 = combinations[,2]
    )
  
  delta_mu <- as.data.frame(rstan::summary(fit,pars="delta_mu")$summary)
  delta_mu <- delta_mu%>%na.omit()
  delta_mu <- cbind(delta_mu_data,parameter="delta_mu",delta_mu[,c("mean","2.5%","97.5%")])
  
  }
  if(t==1){
    delta_mu <- "only possible for more than one time point"
  }
    
  # combine results in one list
  result <- list(mu = mu, 
                 sigma = sigma , 
                 lambda = lambda, 
                 delta_mu= delta_mu, 
                 euclidean_distance=distances)

  # Map metabolite IDs to KEGG and names
  metabolites <- unique(data_df[c(metabolite, kegg)])

  # # Join metadata: for all list elements
  result <- lapply(result,function(x){
    if(is.data.frame(x)&metabolite%in%colnames(x)){
      x <- left_join(x,metabolites,by=metabolite)
      x <- x %>% relocate (!!kegg, .after=!!metabolite)
    }
    return(x)
  })
  # result <- left_join(result, metabolites, by = "metabolite")
  # result <- result %>% relocate("metabolite", .after = "metabolite")
  # result <- result %>% relocate(kegg, .after = "metabolite")
  # result <- unique(result)

  # if input is a SummarizedExperiment object, store estimates in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["estimates_dynamics"]] <- result
    return(data)
  } else {
    # otherwise, return the list of fits
    return(result)
  }
}

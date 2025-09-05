#' cluster dynamics profiles of metabolites
#'
#' Convenient wrapper function for clustering of metabolite dynamics
#' employing the "hybrid" method of the \link[dynamicTreeCut]{dynamicTreeCut}
#' package for clustering and \link[stats]{hclust} for computing of distance
#' matrix and hierarchical clustering needed as input for dynamicTreeCut.
#' Provides bootstrapping of clustering solution from posterior estimates of the model.
#'
#' @param data data frame or colData of a \link[SummarizedExperiment]{SummarizedExperiment}  used used to fit dynamics model
#' @param fit model fit obtained by fit_dynamics_model(). Needed if data is not a SummarizedExperiment object for which the model fit is saved in metadata[["dynamic_fit"]]
#' @param estimates output of estimates_dynamics function, needed if data is not a SummarizedExperiment object for which the model estimates are saved in metadata[["estimates_dynamics"]]
#' @param distance distance method to be used as input for hierarchical clustering \link[stats]{dist}
#' can be "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#' @param agglomeration agglomerative method to be used for hierarchical clustering \link[stats]{hclust}
#' can be "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid"
#' @param minClusterSize minimum number of metabolites per of cluster \link[dynamicTreeCut]{cutreeDynamic}
#' @param deepSplit rough control over sensitivity of cluster analysis. Possible values are 0:4,
#' the higher the value, the more and smaller clusters will be produced by \link[dynamicTreeCut]{cutreeDynamic}
#' @param B number of bootstraps
#'
#' @returns a list with dataframes named by experimental condition or
#' if data is a \link[SummarizedExperiment]{SummarizedExperiment} object clustering
#' results are stored in metadata under "cluster"
#'
#' @export
#'
#' @seealso [fit_dynamics_model()], [estimates_dynamics()], [plot_cluster()]
#'
#' @import SummarizedExperiment
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom S4Vectors metadata
#' @importFrom rstan extract
#' @import tidyr
#' @import dplyr
#' @import ape 
#'
#' @examples
#' data("longitudinalMetabolomics")
#' data <- longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A" &
#'   longitudinalMetabolomics$metabolite %in% c("ATP", "L-Alanine", "GDP")]
#' data <- fit_dynamics_model(
#'   data = data,
#'   scaled_measurement = "m_scaled", assay = "scaled_log",
#'   max_treedepth = 14, adapt_delta = 0.95, iter = 2000, cores = 1, chains = 1
#' )
#' data <- estimates_dynamics(
#'   data = data
#' )
#' data <- cluster_dynamics(data)
#' S4Vectors::metadata(data)[["cluster"]][["A"]][["data"]]
#'
cluster_dynamics <- function(data, fit, distance = "euclidean",
                             agglomeration = "ward.D2",
                             minClusterSize = 1, deepSplit = 2,
                             B = 1000) {
 
  # Input checks
  if (!is.list(estimates)&!inherits(data, "SummarizedExperiment")){
    stop("'data' must be a SummarizedExperiment object or provide estimates")
  }
  # check input class and convert SummarizedExperiment to dataframe
  if (is(data, "SummarizedExperiment")) {
    estimates <- metadata(data)[["estimates_dynamics"]]
    fit <- metadata(data)[["dynamic_fit"]]
  }
  if (!is.list(estimates)) {
    stop("'estimates' must be a list of dataframes obtained by estimates_dynamics()")
  }
  if(!inherits(fit,"stanfit")){
    stop("'fit' must be a modelfit obtained by fit_dynamics_model()")
  }
   
  # Input checks !!!!
  # 
  # check for stanfit
  # check for estimates
  # check for data
  # check for >1 time point per metabolite and condition
  # check for column names metabolite, time, condition
  
  # # check for correct number of mu_mean for every time.ID and metabolite combination
  #   time_ids <- unique(data_df$time)
  #   metabolites <- unique(data$metabolite)
  #   for (metabolite in metabolites) {
  #     if (nrow(unique(data[data$metabolite == metabolite, ])) != length(time_ids)) {
  #       stop("dataframes in 'data' must have a mu_mean value for every time.ID per metabolite'")
  #     }
  #   }
  # }

  # attach variables to function
  metabolite <- NULL
  time.ID <- NULL
  mu_mean <- NULL

  # two steps
  # 1) on mean estimates get clustering solution
  # 2) bootstrapping on clustering solution for bubbletree


# 1) on mean estimates get clustering solution
  # get estimates of mu (scaled metabolite abundance per condition and metabolite)
    mu <- estimates[["mu"]]
    # only select mean estimates
    mu <- mu%>%select(!!metabolite,!!kegg,!!time,!!condition,mean)%>%
      mutate(time=paste0(time,"_mean"))%>% # add label to mean 
      pivot_wider(names_from = !!time, values_from = mean)
    
    # split per condition
    mu_split <- split.data.frame(mu,f=as.factor(mu[[condition]]))
    # get mean clustering solution
    cluster_mean <- lapply(X = mu_split,FUN = .hierarchical_clustering,
                           distance = distance,
                           agglomeration = agglomeration,
                           minClusterSize = minClusterSize,
                           deepSplit = deepSplit)

    
    # bootstrapping with phylograms (function adapted from snaketron/cellmig)
    # get posterior and make different conditions accesible
    posterior <- as.data.frame(rstan::extract(fit,pars="mu"))
    names <- names(posterior)
    x <- as.data.frame(do.call(rbind,strsplit(x=names,split="[.]")))
    colnames(x) <- c("parameter",metabolite,time,condition)
    mu <- estimates[["mu"]]
    data_names <- mu%>%select(!!metabolite,!!time,!!condition)%>%distinct()
    new_names <- paste0("mu.",data_names[[metabolite]],".",data_names[[time]],".",
                        data_names[[condition]])
    names(posterior) <- new_names
    
    # turn into useable format for bootstrapping
    draws <- nrow(posterior)
    posterior$draw <- seq_len(draws)
    posterior <- posterior%>%pivot_longer(cols=-draw,names_to = "ID", values_to = "posterior")
    x <- do.call(rbind,strsplit(posterior$ID,"[.]"))
    posterior$parameter <- x[,1]
    posterior[[metabolite]] <- x[,2]
    posterior[[time]] <- x[,3]
    posterior[[condition]] <- x[,4]
    posterior <- posterior%>%select(-ID)%>%
      pivot_wider(names_from = !!time,values_from = posterior)
    posterior_split <- split.data.frame(posterior,posterior[[condition]])
    
    # get bootstrapping of clustering solution
    boot_ph <- lapply(X = posterior_split,
                      FUN = .get_boot_ph,
                      distance = distance, 
                      agglomeration = agglomeration, 
                      B = B)
    # get clades 
    clades <- list()
    for (i in unique(mu[[condition]])){
      clades[[i]] <- prop.clades(phy = cluster_mean[[i]]$mean_phylo, 
                                 x = c(boot_ph[[i]]), part = NULL,
                  rooted = is.rooted(cluster_mean[[i]]$mean_phylo))
      cluster_mean[[i]]$mean_phylo$nodel.label <- clades[[i]]
      if(all(cluster_mean[[i]]$mean_phylo$tip.label==as.numeric(as.factor(cluster_mean[[i]]$data[[metabolite]])))){
        cluster_mean[[i]]$mean_phylo$tip.label <- as.numeric(as.factor(cluster_mean[[i]]$data[[metabolite]]))
      }
      na_nodes <- which(is.na(cluster_mean[[i]]$mean_phylo$nodel.label))
      if(length(na_nodes)!=0){
        cluster_mean[[i]]$mean_phylo$nodel.label[na_nodes] <- 0
      }
    }

    
    # if(all(main_ph$tip.label == !!metabolite!!!)) {
    #   mean_ph$tip.label <- gns
    # }
    
    # b = 0 for these nodes
    na_nodes <- which(is.na(mean_ph$node.label))
    if(length(na_nodes)!=0) {
      mean_ph$node.label[na_nodes] <- 0
    }
    
    # output
    # 1) dataframe with mean clustering solution -> ORA/compare functions
    # 2) plot_cluster: bubbletree data, dataframe with mean clustering solution -> lineplots
    # 3) vector/dataframe of tip order (needed for ORA to align?)

  
    
    result <- list(cluster_mean=cluster_mean,cluster_bootstrap=cluster_bootstrap)
  # if input is a SummarizedExperiment object, store the fits in the metadata
  if (is(data, "SummarizedExperiment")) {
    metadata(data)[["cluster"]] <- result
    return(data)
  } else {
    # otherwise, return the list of fits
    return(result)
  }
}

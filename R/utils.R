#' input checks for fit_dynamics_model
#' @keywords internal
.check_fit_dynamics_input <- function(model, data, metabolite,
                                     time, condition,
                                     scaled_measurement,
                                     counts, assay, chains, cores, adapt_delta,
                                     max_treedepth, iter, warmup){
  
  if(!model%in%c("scaled_log","raw_plus_counts")){
    stop("'model' must be either 'scaled_log' or 'raw_plus_counts'")
  }
  if(model=="scaled_log"){
    # hint user if data is standardized
    message("Are your metabolite concentrations normalized and standardized?
          We recommend normalization by log-transformation.
          Scaling and centering (mean=0, sd=1) should be metabolite and condition specific.")
  }
  if(model=="raw_plus_counts"){
    # hint user if data is standardized
    message("Are your metabolite concentrations non-normalized and non-scaled?
          The chosen 'raw_plus_counts' model expects raw metabolite concentrations.")
    
    # Input checks
    if (!is.data.frame(counts)) {
      stop("'counts' must be a dataframe if you chose model 'raw_plus_counts'.
                If you cannot provide cell counts choose model 'scaled_log'")
    }
    if (!all(c(time, condition, "counts") %in% colnames(counts))) {
      stop("'counts' must contain columns named 'time','condition', and 'counts'")
    }
  }
  
  # Input checks
  if (!is.data.frame(data) && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a dataframe or a SummarizedExperiment object")
  }
  # check if all input variables are positive integers
  if (!all(vapply(c(iter, warmup, max_treedepth), function(x) is.numeric(x) && x > 0 && x %% 1 == 0, logical(1)))) {
    stop("'iter', 'warmup', and 'max_treedepth' must be positive integers")
  }
  if (!is.numeric(adapt_delta) | !(adapt_delta > 0 & adapt_delta < 1)) {
    stop("'adapt_delta' must be numeric and in the range (0;1)")
  }
  # check if all input variables are character vectors
  if (!all(vapply(list(metabolite, time, condition, scaled_measurement), is.character, logical(1)))) {
    stop("'metabolite', 'time', 'condition', and 'scaled_measurement' must be a character vector specifying a column name of data")
  }
}

#' euclidean distance
#' compare_dynamics()
#' @param a a numeric vector
#' @param b a numeric vector of same length as a
#' @importFrom stats dist
#' @return euclidean distance between vectors
#' @keywords internal

.eu <- function(a, b) {
  temp <- rbind(a, b)
  dist <- stats::dist(temp, method = "euclidean")
  return(dist)
}

#' Jaccard Index: intersection/union
#' compare_metabolites()
#' @param a a vector
#' @param b a vector
#'
#' @return Jaccard Index of a and b
#' @keywords internal
#'
.similarity <- function(a, b) {
  # test which vector is bigger and compare smaller to bigger
  if (length(a) < length(b)) {
    # intersection
    temp <- is.element(a, b) == TRUE
  } else {
    temp <- is.element(b, a) == TRUE
  }
  intersection <- length(temp[temp == TRUE])
  # union=unique metabolites per set + intersection
  sim <- intersection / sum(
    (length(a) - intersection),
    (length(b) - intersection), intersection
  )
  return(sim)
}


#' get bootstrapps for clustering of dynamics vectors (cluster_dynamics function)
#'
#' @keywords internal
.get_boot_tp <- function(x, distance, method, B) {
  
  # !!! all for mean form summary(fit,"parameter")$summary
  
  # eg <- x$posteriors$delta_t  ### adapt to our data set
  # 
  # if(missing(groups)==FALSE) {
  #   if(any(!groups %in% unique(eg$group))) {
  #     stop("selected treatment groups not found in data")
  #   }
  #   eg <- eg[eg$group %in% groups, ]
  # }
  # 
  # # hclust -> main tree
  eg <- x
  gs <- as.numeric(as.factor(eg$metabolite))
  gns <- eg$metabolite
  v <- eg%>%group_by(metabolite,radiation.dose,cell.line)%>%
    mutate(mean_1h=mean(`1h`),mean_3h=mean(`3h`),mean_6h=mean(`6h`),mean_24h=mean(`24h`))%>%
    select(-c(`1h`:`24h`,draw:chain))%>%distinct()
  v <- as.data.frame(v)
  colnames(v) <- gsub("mean_","",colnames(v))
  rownames(v) <- v$metabolite
  # names(v) <- eg$group_id
  hc <- hclust(dist(as.matrix(v[,-c(1:4)]), method = distance), method = method)
  main_ph <- as.phylo(x = hc)
  
  
  # extract posterior
  # e <- extract(x$f, par = "mu_group")$mu_group[, gs]  # extracts 
  
  e <- x
  # get B samples from posterior
  samples <- sample(x= seq_len(draws),size=min(nrow(e),B),replace = TRUE)
  
  boot_ph <- c()
  for(i in seq_len(length(samples))) {
    
    # hclust
    temp <- e[which(e$draw==samples[i]),c(5,8:11)] # metabolite + time vector
    rownames(temp) <- temp$metabolite
    hc <- hclust(dist(as.matrix(temp[,-1]), method = distance), method = method)
    ph <- as.phylo(x = hc)
    
    if(i == 1) {
      boot_ph <- ph
    }
    else {
      boot_ph <- c(boot_ph, ph)
    }
  }
}



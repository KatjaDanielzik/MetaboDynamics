check_fit_dynamics_input <- function(model, data, metabolite,
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

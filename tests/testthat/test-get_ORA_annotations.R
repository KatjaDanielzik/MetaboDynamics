test_that("get_ORA_annotations: input checks", {
  
  # Mock data for testing
  data_df <- as.data.frame(cbind(
    KEGG = c("C00031", "C00123", "C00267"),
    metabolite = c("Glucose", "Fructose", "Sucrose")
  ))
  
  # SummarizedExperiment example
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:6, ncol = 3)),
    colData = data_df
  )
  
  # Check 1: update_background must be a logical (boolean) value
  expect_error(get_ORA_annotations(data = data_df, update_background = "yes"),
               "'update_background' must be either 'TRUE' or 'FALSE'")

  # Check 2: data must be a dataframe or SummarizedExperiment object
  expect_error(get_ORA_annotations(data = matrix(1:9, ncol = 3)), 
               "'data' must be a dataframe or colData of a SummarizedExperiment object")
  
  # Check 3: kegg and metabolite_name must be character vectors
  expect_error(get_ORA_annotations(data = data_df, kegg = "KEGG", metabolite_name = 1), 
               "'metabolite_name' must be a character vector specifying a column name of data")
  expect_error(get_ORA_annotations(data = data_df, kegg = 1, metabolite_name="metabolite"), 
               "'kegg' must be a character vector specifying a column name of data")

  
  # Check 4: data must contain columns named 'kegg' and 'metabolite_name'
  expect_error(get_ORA_annotations(data_df, kegg = "nonexistent", metabolite_name = "metabolite"), 
               "'data' must contain columns named 'kegg', and 'metabolite_name'")
  expect_error(get_ORA_annotations(data_df, kegg = "KEGG", metabolite_name = "nonexistent"), 
               "'data' must contain columns named 'kegg', and 'metabolite_name'")
  
})
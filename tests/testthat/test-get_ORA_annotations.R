# Mock data for testing
mock_data <- as.data.frame(cbind(
  KEGG = c("C00041", "C00362", "C00360"),
  metabolite = c("L-Alanine", "dGMP", "dAMP")
))

test_that("get_ORA_annotations: input checks", {
  # SummarizedExperiment example
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:6, ncol = 3)),
    colData = mock_data
  )

  # update_background must be a logical (boolean) value
  expect_error(
    get_ORA_annotations(data = mock_data, update_background = "yes"),
    "'update_background' must be either 'TRUE' or 'FALSE'"
  )

  # data must be a dataframe or SummarizedExperiment object
  expect_error(
    get_ORA_annotations(data = matrix(1:9, ncol = 3)),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # kegg and metabolite_name must be character vectors
  expect_error(
    get_ORA_annotations(data = mock_data, kegg = "KEGG", metabolite_name = 1),
    "'metabolite_name' must be a character vector specifying a column name of data"
  )
  expect_error(
    get_ORA_annotations(data = mock_data, kegg = 1, metabolite_name = "metabolite"),
    "'kegg' must be a character vector specifying a column name of data"
  )

  #  data must contain columns named 'kegg' and 'metabolite_name'
  expect_error(
    get_ORA_annotations(mock_data, kegg = "nonexistent", metabolite_name = "metabolite"),
    "'data' must contain columns named 'kegg', and 'metabolite_name'"
  )
  expect_error(
    get_ORA_annotations(mock_data, kegg = "KEGG", metabolite_name = "nonexistent"),
    "'data' must contain columns named 'kegg', and 'metabolite_name'"
  )
})


test_that("get_ORA_annotations:output_checks", {
  result <- get_ORA_annotations(data = mock_data, kegg = "KEGG", metabolite_name = "metabolite", update_background = FALSE)
  annotation_df <- result[["annotation"]]

  expect_true("KEGG" %in% colnames(annotation_df))
  expect_true("module_id" %in% colnames(annotation_df))
  expect_true("module_name" %in% colnames(annotation_df))

  # Check that the function handles missing KEGG entries gracefully
  mock_data_with_na <- data.frame(
    KEGG = c("C00041", NA, "C00360"),
    metabolite = c("ATP", "ADP", "AMP")
  )
  result <- get_ORA_annotations(data = mock_data_with_na, kegg = "KEGG", metabolite_name = "metabolite", update_background = FALSE)
  annotation_df <- result[["annotation"]]

  # Expect the annotation dataframe to handle NA gracefully, without errors
  expect_equal(nrow(annotation_df), 2)
})

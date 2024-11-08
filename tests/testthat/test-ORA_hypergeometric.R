test_that("ORA_hypergeometric: input checks", {
  
  # Create valid mock data for testing
  valid_background <- data.frame(KEGG = c("k1", "k2"), middle_hierarchy = c("mh1", "mh2"))
  valid_annotations <- data.frame(KEGG = c("k1", "k2"), middle_hierarchy = c("mh1", "mh2"))
  valid_clusters <- data.frame(KEGG = c("k1", "k2"), cluster = c("A", "B"))
  valid_tested_column <- "middle_hierarchy"
  
  # Invalid 'background' input (non-dataframe input)
  invalid_background <- list(KEGG = c("k1", "k2"))
  expect_error(ORA_hypergeometric(background = invalid_background, 
                                  annotations = valid_annotations,
                                  clusters <- valid_clusters,
                                  tested_column = valid_tested_column),
   "'background' must be a dataframe obtained by get_ORA_annotations()")
  
  # Invalid 'annotations' input (non-dataframe input)
  invalid_annotations <- list(KEGG = c("k1", "k2"))
  expect_error(ORA_hypergeometric(background = valid_background, 
                                  annotations = invalid_annotations,
                                  clusters <- valid_clusters,
                                  tested_column = valid_tested_column),
           "'annotations' must be a dataframe obtained by get_ORA_annotations()")

  # Invalid 'clusters' input (non-dataframe input)
  invalid_clusters <- list(KEGG = c("k1", "k2"), cluster = c("A", "B"))
  expect_error(ORA_hypergeometric(background = valid_background, 
                                  annotations = valid_annotations,
                                  clusters <- invalid_clusters,
                                  tested_column = valid_tested_column),
                                  "'clusters' must be a dataframe")
  
  #  Missing required columns in 'clusters'
  missing_columns_clusters <- data.frame(KEGG = c("k1", "k2"))
  expect_error(ORA_hypergeometric(background = valid_background, 
                                  annotations = valid_annotations,
                                  clusters <- missing_columns_clusters,
                                  tested_column = valid_tested_column), 
               "'clusters' must contains columns 'KEGG' and 'cluster'")

  # Invalid 'tested_column' input (not a character)
  invalid_tested_column <- 123
  expect_error(ORA_hypergeometric(background = valid_background, 
                                  annotations = valid_annotations,
                                  clusters <- valid_clusters,
                                  tested_column = invalid_tested_column), 
               "'tested_column' must be a character vector")
  
  # 'tested_column' not present in 'clusters', 'background', or 'annotations'
  invalid_tested_column <- "invalid_column"
  expect_error(ORA_hypergeometric(background = valid_background, 
                                  annotations = valid_annotations,
                                  clusters <- valid_clusters,
                                  tested_column = invalid_tested_column), 
               "'tested_column' must be a column of 'clusters', 'background' and 'annotations'")

  # Partially missing 'tested_column' in 'clusters', 'background', or 'annotations'
  missing_column_background <- data.frame(KEGG = c("k1", "k2"))
  expect_error(ORA_hypergeometric(background = missing_column_background, 
                                  annotations = valid_annotations,
                                  clusters <- valid_clusters,
                                  tested_column = invalid_tested_column), 
               "'tested_column' must be a column of 'clusters', 'background' and 'annotations'")
})
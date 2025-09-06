# Create valid mock data for testing
valid_annotations <- data.frame(
  KEGG = c("C00001", "C00002"),
  module_id = c(1, 2),
  module_name = c("Module1", "Module2"),
  middle_hierarchy = c("MH1", "MH2")
)
valid_background <- data.frame(
  KEGG = c("C00001", "C00002"),
  module_id = c(1, 2),
  module_name = c("Module1", "Module2"),
  middle_hierarchy = c("MH1", "MH2")
)
valid_clusters <- list(A=list(data=data.frame(
  KEGG = c("C00001", "C00002"),
  cluster = c(1, 2),
  condition = c("A", "B"))
))

invalid_background <- data.frame(
  KEGG = c("C00001", "C00002"),
  module_id = c(1, 2),
  module_name = c("Module1", "Module2")
)

valid_tested_column <- "middle_hierarchy"

test_that("ORA_hypergeometric: input checks", {
  # Invalid 'background' input (non-dataframe input)
  invalid_background <- list(KEGG = c("k1", "k2"))
  expect_error(
    ORA_hypergeometric(
      background = invalid_background,
      annotations = valid_annotations,
      data <- valid_clusters,
      tested_column = valid_tested_column
    )
  )

  # Invalid 'annotations' input (non-dataframe input)
  invalid_annotations <- list(KEGG = c("k1", "k2"))
  expect_error(
    ORA_hypergeometric(
      background = valid_background,
      annotations = invalid_annotations,
      data <- valid_clusters,
      tested_column = valid_tested_column
    )
  )

  #  Missing required columns in 'data'
  missing_columns_clusters <- list(A=list(data=data.frame(
    data = c(1,2),
    cluster = c(1, 2),
    condition = c("A", "B"))
  ))
  expect_error(
    ORA_hypergeometric(
      background = valid_background,
      annotations = valid_annotations,
      data = missing_columns_clusters,
      tested_column = valid_tested_column
    ),
    "'data' must contains columns 'KEGG', 'cluster' and 'condition'"
  )

  # Invalid 'tested_column' input (not a character)
  expect_error(
    ORA_hypergeometric(
      background = valid_background,
      annotations = valid_annotations,
      data <- valid_clusters,
      tested_column = 123
    ),
    "'tested_column' must be a character vector"
  )

  # 'tested_column' not present in 'data', 'background', or 'annotations'
  invalid_tested_column <- "invalid_column"
  expect_error(ORA_hypergeometric(
    background = valid_background,
    annotations = valid_annotations,
    data <- valid_clusters,
    tested_column = invalid_tested_column
  ))
})


test_that("ORA_hypergeometric:output_checks", {
  # Ensure that the results contain the OvE values and related columns
  result <- ORA_hypergeometric(
    background = valid_background,
    annotations = valid_annotations,
    data = valid_clusters,
    tested_column = "middle_hierarchy"
  )
  expect_true("condition" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))
  expect_true("OvE_gen" %in% colnames(result))
  expect_true("OvE_gen_lower" %in% colnames(result))
  expect_true("OvE_gen_higher" %in% colnames(result))
  expect_true("OvE_gen_median" %in% colnames(result))
  # Example check: ensure that OvE_gen is a numeric value
  expect_true(is.numeric(result$OvE_gen))
  # Ensure OvE_gen values are meaningful, non-zero for observed over expected.
  expect_true(all(result$OvE_gen > 0))
})

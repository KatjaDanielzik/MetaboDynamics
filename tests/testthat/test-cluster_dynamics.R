test_that("cluster_dynamics: input checks", {
  # Function returns an error when data is not a list of dataframes or SummarizedExperiment object
  test_data <- data.frame(a = c(1, 2, 3))
  expect_error(cluster_dynamics(test_data), "'data' must be a list of dataframes or a SummarizedExperiment object")

  # Function returns an error when data is a list and contains non-dataframe elements
  test_data <- list(
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = c(1, 2, 3),
      time.ID = c(1, 2, 3)
    ),
    matrix(1:9, nrow = 3, ncol = 3)
  )
  expect_error(cluster_dynamics(test_data), "'data' must be a list of dataframes if it is not a SummarizedExperiment object")

  # Function returns an error when dataframes in data do not contain required columns
  test_data <- list(
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      time.ID = c(1, 2, 3)
    ),
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = c(4, 5, 6),
      time.ID = c(1, 2, 3)
    )
  )
  expect_error(
    cluster_dynamics(test_data),
    "dataframes in 'data' must contain columns named 'metabolite', 'mu_mean' and 'time.ID'"
  )

  # Function returns an error when there are missing mu_mean values for
  # any metabolite-time.ID combination
  test_data <- list(
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = c(1, 2, 3),
      time.ID = c(1, 2, 3)
    ),
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = c(4, 5, 6),
      time.ID = c(1, 2, 3)
    )
  )
  expect_error(
    cluster_dynamics(test_data),
    "dataframes in 'data' must have a mu_mean value for every time.ID per metabolite'"
  )
})

test_that("cluster_dynamics: output checks", {
  # Function returns a list with the correct number of elements when data is
  # a list of data frames
  test_data <- list(
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = rep(c(1, 2, 3), 3),
      time.ID = rep(c(1, 2, 3), each = 3)
    ),
    data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      mu_mean = rep(c(1, 2, 3), 3),
      time.ID = rep(c(1, 2, 3), each = 3)
    )
  )
  result <- cluster_dynamics(test_data)
  expect_equal(length(result), 2)

  # Function returns correct named list elements and elements have correct
  # classes    expect_true(all(c("data", "distance_matrix", "tree") %in% names(result)))
  result <- result[[1]]
  expect_true(is.data.frame(result$data))
  expect_true(is.matrix(result$distance_matrix))
  expect_true(inherits(result$tree, "hclust"))
})

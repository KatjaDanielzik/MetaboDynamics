test_that("plot_PPC: input checks", {
  
  # Mock data for testing
  valid_data <- data.frame(
    time = c(0, 1, 2, 3),
    condition = c("A", "B", "A", "B"),
    scaled_measurement = c(0.5, 1.2, 0.8, 1.5)
  )
  
  valid_posterior <- data.frame(
    mu_mean = c(0.1, 0.2),
    mu_lower = c(0.05, 0.1),
    mu_higher = c(0.15, 0.3)
  )
  
  invalid_data <- list(
    time = c(0, 1, 2, 3),
    condition = c("A", "B", "A", "B")
  )
  
  invalid_posterior <- list(
    mu_mean = c(0.1, 0.2),
    mu_lower = c(0.05, 0.1)
  )
  
  # Check if 'data' is a dataframe or SummarizedExperiment
  expect_error(
    plot_PPC(posterior = valid_posterior, data = invalid_data, 
             scaled_measurement = "scaled_measurement"),
    "'data' must be a dataframe or colData of a SummarizedExperiment object"
  )
  
  # Check if 'posterior' is a dataframe
  expect_error(
    plot_PPC(invalid_posterior, valid_data, "scaled_measurement"),
    "'posterior' must be a dataframe obtained by diagnostics_dynamics()"
  )
  
  # Check if 'scaled_measurement' is a character vector
  expect_error(
    plot_PPC(valid_posterior, valid_data, 123),
    "'scaled_measurement' must be a character vector specifying a column name of data"
  )
  
  # Check if 'scaled_measurement' exists in 'data' column names
  expect_error(
    plot_PPC(valid_posterior, valid_data, "non_existing_column"),
    "'data' must contain a column named 'scaled_measurement'")
  })
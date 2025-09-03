test_that("plot_PPC: input checks", {
  # Mock data for testing
  valid_data <- data.frame(
    time = c(0, 1, 2, 3),
    condition = c("A", "B", "A", "B"),
    scaled_measurement = c(0.5, 1.2, 0.8, 1.5),
    metabolite = rep("A", 4)
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
    plot_PPC(
      posterior = valid_posterior, data = invalid_data,
      scaled_measurement = "scaled_measurement"
    ),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # Check if 'posterior' is a list
  expect_error(
    plot_PPC(
      posterior = invalid_posterior, data = valid_data,
      "scaled_measurement"
    ),
    "'posterior' must be a data frame obtained by diagnostics_dynamics()"
  )

  # Check if 'scaled_measurement' is a character vector
  expect_error(
    plot_PPC(valid_posterior, valid_data, scaled_measurement = 123),
    "'scaled_measurement' must be a character vector specifying a column name of data"
  )

  # Check if 'scaled_measurement' exists in 'data' column names
  expect_error(
    plot_PPC(valid_posterior, valid_data, scaled_measurement = "non_existing_column"),
    "'data' must contain a column named 'scaled_measurement' and 'metabolite'"
  )
})

test_that("plot_PPC:output_checks", {
  posterior_example <- data.frame(
    time.ID = rep(1:5, each = 10),
    condition = "A",
    posterior = rnorm(50)
  )

  # Example data dataframe (from SummarizedExperiment)
  data_example <- data.frame(
    metabolite = rep("ATP", 10),
    time = rep(1:5, each = 2),
    m_scaled = rnorm(10)
  )
  plot <- plot_PPC(
    posterior = posterior_example, data = data_example,
    scaled_measurement = "m_scaled"
  )
  expect_s3_class(plot, "gg")
})

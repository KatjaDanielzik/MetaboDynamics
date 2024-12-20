test_that("plot_estimates: input checks", {
  # Mock valid data for testing
  valid_data <- data.frame(time = c(1, 2, 3), metabolite = c("A", "B", "C"))
  valid_estimates <- data.frame(
    condition = c("A"),
    metabolite.ID = c(1, 2),
    metabolite = c("A", "B"),
    higher = c(1.5, 2.5),
    lower = c(0.5, 1.5)
  )
  valid_estimates <- list(valid_estimates)

  # Invalid data mocks
  invalid_data <- list(time = c(1, 2, 3))
  invalid_estimates <- as.data.frame(cbind(condition = c("dose1", "dose2")))

  # 'data' must be a dataframe or SummarizedExperiment object
  expect_error(
    plot_estimates(estimates = valid_estimates, data = invalid_data),
    "'data' must be a dataframe or colData of a SummarizedExperiment object"
  )

  # delta_t' must be logical
  expect_error(
    plot_estimates(estimates = valid_estimates, data = valid_data, delta_t = "yes"),
    "'delta_t' must be either 'TRUE' or 'FALSE'"
  )

  # 'dynamics' must be logical (TRUE or FALSE)
  expect_error(
    plot_estimates(estimates = valid_estimates, data = valid_data, dynamics = 1),
    "'dynamics' must be either 'TRUE' or 'FALSE'"
  )
})

test_that("plot_estimates:output_checks", {
  # Example data
  A <- data.frame(
    condition = rep("A", 4),
    metabolite.ID = 1,
    metabolite = rep("ATP", 4),
    time.ID = 1:4,
    mu_mean = rnorm(4),
    mu_lower = rnorm(4),
    mu_higher = rnorm(4),
    delta_mu_mean = c(rnorm(3), NA),
    delta_mu_lower = c(rnorm(3), NA),
    delta_mu_higher = c(rnorm(3), NA)
  )
  estimates <- list(A = A, B = A)

  data <- data.frame(
    time = rep(1:2, each = 10),
    condition = rep(c("A", "B"), each = 10),
    metabolite = rep("ATP", 20)
  )


  plots <- plot_estimates(estimates = estimates, data = data, delta_t = TRUE, dynamics = TRUE)

  expect_true("plot_timepoint_differences" %in% names(plots)) # Ensure timepoint differences plot is still generated
  expect_true("plot_dynamics" %in% names(plots))
})

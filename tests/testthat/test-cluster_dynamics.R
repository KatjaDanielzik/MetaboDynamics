# Define a dummy data frame
dummy_data <- data.frame(
  time = rep(1:2, each = 2 * 3),
  metabolite = rep(letters[1:2], times = 2 * 3),
  condition = rep("A", length.out = 12),
  scaled_measurement = rnorm(12),
  kegg = rep(rnorm(2), times = 2 * 3)
)

# Define a dummy fit object
dummy_fit <- fit_dynamics_model(model = "scaled_log",data=dummy_data,
                                scaled_measurement="scaled_measurement",
                                cores=1,
                                chains=1)

# Define a dummy estimates list
dummy_estimates <- estimates_dynamics(data=dummy_data,fit=dummy_fit,kegg="kegg")

# Define the tests
test_that("cluster_dynamics:input_checks", {
  # Test that the function throws an error if data is not a SummarizedExperiment object or if estimates are not provided
  expect_error(
    cluster_dynamics(data = "not a SummarizedExperiment object", fit = dummy_fit),
    "'data' must be a SummarizedExperiment object or provide estimates"
  )
  
  # Test that the function throws an error if estimates is not a list of data frames obtained by estimates_dynamics()
  expect_error(
    cluster_dynamics(data = dummy_data, fit = dummy_fit, estimates = "not a list of data frames"),
    "'data' must be a SummarizedExperiment object or provide estimates"
  )
  
  # Test that the function throws an error if mu is not a data frame
  expect_error(
    cluster_dynamics(data = dummy_data, fit = dummy_fit, estimates = list(mu = "not a data frame")),
    "'mu' must be a data frame"
  )
  
  # Test that the function throws an error if fit is not a modelfit obtained by fit_dynamics_model()
  expect_error(
    cluster_dynamics(data = dummy_data, fit = "not a modelfit", estimates = dummy_estimates),
    "'fit' must be a modelfit obtained by fit_dynamics_model()"
  )
  
  # Test that the function throws an error if every metabolite and condition does not have at least one time point for clustering
  test <- dummy_estimates
  test$mu$time[2] <- NA
  expect_error(
    cluster_dynamics(data = dummy_data, fit = dummy_fit, estimates = test),
    "every metabolite and condition must have at least one time point for clustering"
  )
  
  # Test that the function throws an error if the number of time points is not the same for all metabolites and conditions
  test <- dummy_estimates
  test$mu$time[2] <- 1
  expect_error(
    cluster_dynamics(data = dummy_data, fit = dummy_fit, estimates = test),
    "number of time points must be the same for all metabolites and conditions"
  )
})

test_that("cluster_dynamics:output_checks", {
  # Test that the function returns a list with the expected elements
  results <- cluster_dynamics(data = dummy_data, fit = dummy_fit, estimates = dummy_estimates)
  expect_type(results, "list")
  expected_elements <- c("cluster_mean", "cluster_bootstrap")
  expect_true(all(expected_elements %in% names(results)))
  
})
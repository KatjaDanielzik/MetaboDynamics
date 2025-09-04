
# Define a dummy data frame
dummy_data <- data.frame(
  time = rep(1:2, each = 2 * 3),
  metabolite = rep(letters[1:2], times = 2 * 3),
  condition = rep("A", length.out = 12),
  scaled_measurement = rnorm(12),
  kegg = rep(rnorm(2), times = 2 * 3)
)

# Define a dummy estimates list
dummy_estimates <- list(
  delta_mu = as.data.frame(cbind(
    metabolite = rep(letters[1:2], each = 3),
    condition = rep("A", length.out = 6),
    timepoint_1 = rep(1:2, each = 3),
    timepoint_2 = rep(2:1, each = 3),
    mean = rnorm(6),
    `2.5%` = rnorm(6),
    `97.5%` = rnorm(6)
  )),
  euclidean_distances = as.data.frame(cbind(
    metabolite = rep(letters[1:2], each = 3),
    condition_1 = rep("A", length.out = 6),
    condition_2 = rep("B", length.out = 6),
    mean = rnorm(6),
    `2.5%` = rnorm(6),
    `97.5%` = rnorm(6)
  )),
  mu = as.data.frame(cbind(
    metabolite = rep(letters[1:2], each = 3),
    condition = rep("A", length.out = 6),
    time = rep(1:2, each = 3),
    mean = rnorm(6),
    `2.5%` = rnorm(6),
    `97.5%` = rnorm(6)
  ))
)

# Define the tests
test_that("plot_estimates:input_checks", {
  # Test that the function throws an error if data is not a data frame or SummarizedExperiment object
  expect_no_error(
    plot_estimates(data = "not a data frame", estimates = dummy_estimates)
  )
  
  # Test that the function throws an error if estimates is not a list of data frames obtained by estimates_dynamics()
  expect_error(
    plot_estimates(data = dummy_data, estimates = "not a list of data frames"),
    "'data' must be a SummarizedExperiment object or provide estimates"
  )
  
  # Test that the function throws an error if delta_t is not a logical value
  expect_error(
    plot_estimates(data = dummy_data, estimates = dummy_estimates, delta_t = "not a logical value"),
    "'delta_t' must be either 'TRUE' or 'FALSE'"
  )
  
  # Test that the function throws an error if dynamics is not a logical value
  expect_error(
    plot_estimates(data = dummy_data, estimates = dummy_estimates, dynamics = "not a logical value"),
    "'dynamics' must be either 'TRUE' or 'FALSE'"
  )
  
  # Test that the function throws an error if distance_conditions is not a logical value
  expect_error(
    plot_estimates(data = dummy_data, estimates = dummy_estimates, distance_conditions = "not a logical value"),
    "'distance_conditions' must be either 'TRUE' or 'FALSE'"
  )
})

test_that("plot_estimates:output_checks", {
  # Test that the function returns a list of ggplot objects
  results <- plot_estimates(data = dummy_data, estimates = dummy_estimates)
  expect_type(results, "list")

  # Test that the function returns a list with the expected elements
  expected_elements <- c("delta_t", "distance_conditions", "dynamcis")
  expect_true(all(expected_elements %in% names(results)))
})
  
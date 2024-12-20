dummy_data <- data.frame(
  time = rep(1:2, each = 2 * 3),
  metabolite = rep(letters[1:2], times = 2 * 3),
  condition = rep("A", length.out = 12),
  scaled_measurement = rnorm(12),
  kegg = rep(rnorm(2), times = 2 * 3)
)

dummy_fit <- fit_dynamics_model(
  data = dummy_data, chains = 1, cores = 1,
  iter = 1000,
  metabolite = "metabolite",
  condition = "condition", time = "time",
  scaled_measurement = "scaled_measurement"
)


# Define tests
test_that("estimates_dynamics:input_checks", {
  expect_error(
    estimates_dynamics("not a dataframe", fits = dummy_fit),
    "'data' must be a dataframe or colData of a SummarizedExperiment object"
  )
  expect_error(
    estimates_dynamics(dummy_data, fits = "not a list of stanfit objects"),
    "'fits' must be a list of stanfit objects"
  )
  expect_error(
    estimates_dynamics(dummy_data, kegg = 1, condition = 1, time = "time", fits = dummy_fit),
    "'time', 'kegg' and 'condition' must be a character vector specifying a column name of data"
  )
  invalid_data <- dummy_data[, -3]
  expect_error(
    estimates_dynamics(invalid_data, fits = dummy_fit),
    "'data' must contain columns named 'condition', 'kegg' and 'time'"
  )
})

test_that("estimates_dynamics:output_checks", {
  results <- estimates_dynamics(
    data = dummy_data, time = "time",
    condition = "condition", kegg = "kegg",
    fits = dummy_fit, samples = 2, iter = 1000,
    chains = 1
  )
  expect_type(results, "list")
  expect_named(results, c("A"))
  expect_true(all(c("metabolite.ID", "time.ID", "condition", "mu_mean", "sigma_mean")
  %in% colnames(results[["A"]])))
  expect_equal(unique(results[["A"]]$condition), "A")

  # test for correct output
  results <- estimates_dynamics(dummy_data,
    condition = "condition", kegg = "kegg",
    fits = dummy_fit, samples = 5, iter = 1000,
    chains = 1
  )
  result_A <- results[["A"]]
  sample_columns <- grep("^mu_sample_", colnames(result_A), value = TRUE)
  expect_length(sample_columns, 5)
  expect_true("lambda_mean" %in% colnames(result_A))
  expect_true("delta_mu_mean" %in% colnames(result_A), )
})

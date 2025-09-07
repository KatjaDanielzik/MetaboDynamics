# Create a mock dataset for testing
set.seed(123)
mock_data <- data.frame(
  metabolite = rep(c("Metabolite_A", "Metabolite_B"), each = 10),
  time = rep(seq(1, 10), times = 2),
  condition = rep(c("Condition_A", "Condition_B"), each = 10),
  m_scaled = rnorm(20, mean = 0, sd = 1)
)

test_that("fit_dynamics_model:input_checks", {
  # Test: data must be a dataframe or SummarizedExperiment
  expect_error(
    fit_dynamics_model(data = list(), scaled_measurement = "m_scaled"),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # missing columns in the data
  expect_error(
    fit_dynamics_model(
      data = mock_data[, -1], # remove the metabolite column
      scaled_measurement = "m_scaled"
    ),
    "'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'"
  )

  # adapt_delta must be in range [0;1]
  expect_error(fit_dynamics_model(
    data = mock_data,
    adapt_delta = 1.5, scaled_measurement = "m_scaled"
  ))

  expect_error(
    fit_dynamics_model(
      data = mock_data, scaled_measurement = "m_scaled"
    ),
    "Input must contain at least three replicates per metabolite,
      time point and experimental condition."
  )
  
  # Test: counts must be a dataframe if model is 'raw_plus_counts'
  expect_error(
    fit_dynamics_model(model = "raw_plus_counts", data = mock_data, counts = list()),
    "'counts' must be a dataframe if you chose model 'raw_plus_counts'."
  )
  
  # Test: counts must contain columns named 'time','condition', and 'counts'
  expect_error(
    fit_dynamics_model(model = "raw_plus_counts", data = mock_data, counts = data.frame(time = 1:10)),
    "'counts' must contain columns named 'time','condition', and 'counts'"
  )
  
})

test_that("fit_dynamics_model:output_checks", {
  # create triplicates
  mock_data <- data.frame(
    metabolite = rep(c("Metabolite_A", "Metabolite_B"), each = 10),
    time = rep(seq(1, 10), times = 2),
    condition = rep(c("Condition_A", "Condition_B"), each = 10),
    m_scaled = rnorm(20, mean = 0, sd = 1)
  )
  mock_data <- rbind(mock_data, mock_data, mock_data)

  # basic function output
  fit <- fit_dynamics_model(
    data = mock_data,
    scaled_measurement = "m_scaled",
    chains = 1,
    cores = 1,
    iter = 100, # Use fewer iterations for testing purposes
    warmup = 20, adapt_delta = 0.8, max_treedepth = 10
  )

  
  # Test: output must be a 'stanfit' object
  expect_true(inherits(fit, "stanfit"))

})


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
  # column names must be character strings
  expect_error(
    fit_dynamics_model(
      data = mock_data,
      metabolite = 123, time = 1, condition = TRUE, scaled_measurement = "m_scaled"
    ),
    "'metabolite', 'time', 'condition', and 'scaled_measurement' must be a character vector specifying a column name of data"
  )
  # missing columns in the data
  expect_error(
    fit_dynamics_model(
      data = mock_data[, -1], # remove the metabolite column
      scaled_measurement = "m_scaled", time = "time", condition = "condition"
    ),
    "'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'"
  )

  expect_error(
    fit_dynamics_model(
      data = mock_data, # remove the metabolite column
      scaled_measurement = "m_scaled", condition = "concentration", time = "time"
    ),
    "'data' must contain columns named 'metabolite','time','condition', and 'scaled_measurement'"
  )

  # adapt_delta must be in range [0;1]
  expect_error(fit_dynamics_model(
    data = mock_data, condition = "condition", time = "time",
    adapt_delta = 1.5, scaled_measurement = "m_scaled"
  ))

  expect_error(
    fit_dynamics_model(
      data = mock_data, scaled_measurement = "m_scaled",
      condition = "condition"
    ),
    "Input must contain at least three measurements per metabolite,
      time and experimental condition."
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
  fits <- fit_dynamics_model(
    data = mock_data,
    metabolite = "metabolite",
    time = "time",
    condition = "condition",
    scaled_measurement = "m_scaled",
    chains = 1,
    cores = 1,
    iter = 100, # Use fewer iterations for testing purposes
    warmup = 20, adapt_delta = 0.8, max_treedepth = 10
  )

  # output must be a list
  expect_type(fits, "list")

  # Check that the list has one fit per condition
  expect_equal(length(fits), length(unique(mock_data$condition)))

  # Check that each element of the list is a 'stanfit' object
  expect_true(all(sapply(fits, function(x) inherits(x, "stanfit"))))

  # Check that the list has correct names (condition names)
  expect_equal(names(fits), unique(mock_data$condition))
})

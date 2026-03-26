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

  # Test: model must be one of the allowed values
  expect_error(fit_dynamics_model(
    data = mock_data,
    model = "invalid_model",
    scaled_measurement = "m_scaled"
  ), "'model' must be either 'scaled_log' or 'raw_plus_counts'")

  # Test: model_option must be one of the allowed values
  expect_error(fit_dynamics_model(
    data = mock_data,
    model_option = "invalid_option",
    scaled_measurement = "m_scaled"
  ), "'model_option' must be either 'sd_per_time_point' or 'sd_per_condition'")

  # Test: model_option = sd_per_time_point requires at least 3 replicates
  mock_data <- data.frame(
    metabolite = rep(c("Metabolite_A", "Metabolite_B"), each = 10),
    time = rep(seq(1, 10), times = 2),
    condition = rep(c("Condition_A", "Condition_B"), each = 10),
    m_scaled = rnorm(20, mean = 0, sd = 1)
  )

  expect_error(
    fit_dynamics_model(
      data = mock_data,
      model = "scaled_log",
      model_option = "sd_per_time_point",
      scaled_measurement = "m_scaled"
    ),
    "Input must contain at least three replicates per metabolite, time point and experimental condition."
  )

  # Test: model_option = sd_per_condition requires at least 2 replicates
  expect_error(
    fit_dynamics_model(
      data = mock_data,
      model = "scaled_log",
      model_option = "sd_per_condition",
      scaled_measurement = "m_scaled"
    ),
    "Input must contain at least two replicates per metabolite,
      time point and experimental condition. Check diagnostics and PPC carefully
      before using estimates!"
  )

  # Test: counts must be a data frame if model is 'raw_plus_counts'
  expect_error(
    fit_dynamics_model(model = "raw_plus_counts", data = mock_data, counts = list()),
    "'counts' must be a dataframe if you chose model 'raw_plus_counts'."
  )

  # Test: counts must contain columns named 'time','condition', and 'counts'
  expect_error(
    fit_dynamics_model(model = "raw_plus_counts", data = mock_data, counts = data.frame(time = 1:10)),
    "'counts' must contain columns named 'time','condition', and 'counts'"
  )

  # Test: time and condition must match between data and counts
  mock_data <- data.frame(
    metabolite = rep(c("Metabolite_A", "Metabolite_B"), each = 6),
    time = rep(seq(1, 2), each = 3),
    condition = rep(c("Condition_A"), 3),
    m_scaled = rnorm(6, mean = 0, sd = 1)
  )

  time_mismatch <- as.data.frame(cbind(unique(mock_data[, c("condition", "time")]),
    counts = 100
  ))
  time_mismatch$time <- 1

  expect_error(
    fit_dynamics_model(
      model = "raw_plus_counts",
      data = mock_data,
      counts = time_mismatch,
      scaled_measurement = "m_scaled"
    ),
    "data and counts must have the same time points"
  )

  condition_mismatch <- data.frame(
    time = 1:10,
    condition = rep("A", 5),
    counts = 100
  )
  expect_error(
    fit_dynamics_model(
      model = "raw_plus_counts",
      data = mock_data,
      counts = condition_mismatch,
      scaled_measurement = "m_scaled"
    ),
    "data and counts must have the same time points"
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

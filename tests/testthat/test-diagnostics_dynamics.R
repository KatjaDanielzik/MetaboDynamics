dummy_data <- data.frame(
  time = rep(1:2, each = 2 * 3),
  metabolite = rep(letters[1:2], times = 2 * 3),
  condition = rep("A", length.out = 12),
  scaled_measurement = rnorm(12)
)

dummy_fit <- fit_dynamics_model(
  data = dummy_data, chains = 1, cores = 1, iter = 200,
  scaled_measurement = "scaled_measurement"
)

test_that("diagnostics_dynamics: input checks", {
  # Invalid data type for 'data'
  expect_error(
    diagnostics_dynamics(data = list(), fit = dummy_fit),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # Invalid data type for 'fits'
  expect_error(
    diagnostics_dynamics(data = dummy_data, fit = list(123)),
    "'fit' must be a stanfit object"
  )

  # Missing required column in 'data'
  expect_error(
    diagnostics_dynamics(data = dummy_data[, -1], fit = dummy_fit),
    "'data' must contain columns named 'metabolite','time', and 'condition'"
  )
})


test_that("diagnostics_dynamics:output_checks", {
  result <- diagnostics_dynamics(
    data = dummy_data,
    fit = dummy_fit,
    iter = 200,
    warmup = 50,
    chains = 1
  )

  # Check the structure of the output
  expect_true(is.list(result))
  expect_named(result, c("model_diagnostics", "posterior"))

  # Check 'model_diagnostics' dataframe
  diagnostics <- result[["model_diagnostics"]]
  expect_true(is.data.frame(diagnostics))
  expect_true(all(c("metabolite.ID", "condition", "divergences", "treedepth_error") %in% colnames(diagnostics)))

  diagnostics <- result[["model_diagnostics"]]

  # Check that divergences are non-negative integers
  expect_true(all(diagnostics$divergences >= 0))
  expect_true(all(diagnostics$divergences %% 1 == 0))

  # Check Rhat and n_eff values for each timepoint
  rhat_cols <- grep("^rhat_mu", colnames(diagnostics), value = TRUE)
  neff_cols <- grep("^neff_mu", colnames(diagnostics), value = TRUE)

  expect_true(all(diagnostics[, rhat_cols] > 0))
  expect_true(all(diagnostics[, neff_cols] > 0))

  posterior <- result[["posterior"]]

  # Check that posterior is a dataframe
  expect_true(is.data.frame(posterior))

  # Check required columns
  expect_true(all(c("parameter", "posterior", "metabolite.ID", "time.ID") %in% colnames(posterior)))

  # Check that all metabolites and timepoints are present
  unique_metabolites <- unique(dummy_data$metabolite)
  unique_times <- unique(dummy_data$time)

  expect_true(all(unique(posterior$metabolite.ID) %in% seq_along(unique_metabolites)))
  expect_true(all(as.numeric(as.character(posterior$time.ID)) %in% unique_times))
})

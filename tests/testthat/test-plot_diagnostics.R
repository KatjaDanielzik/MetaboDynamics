# Mock valid data for testing
valid_data <- data.frame(
  time = rep(1:2, 10),
  condition = rep(c("A", "B"), each = 10),
  metabolite = rep("ATP", 20)
)
valid_diagnostics <- data.frame(
  metabolite.ID = (1:20), condition = sample(c("A", "B"), 20, replace = TRUE),
  divergences = rnorm(20), treedepth_error = rnorm(20),
  rhat_mu_mean_1 = rnorm(20), rhat_mu_mean_2 = rnorm(20),
  neff_mu_mean_1 = rnorm(20), neff_mu_mean_2 = rnorm(20)
)

test_that("plot_diagnostics: input checks", {
  # Invalid data mocks
  invalid_data <- list(time = c(1, 2, 3))
  invalid_diagnostics <- list(condition = c("dose1", "dose2"))

  # 'data' must be a dataframe or SummarizedExperiment object
  expect_error(
    plot_diagnostics(data = invalid_data, diagnostics = valid_diagnostics),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # 'diagnostics' must be a dataframe
  expect_error(
    plot_diagnostics(data = valid_data, diagnostics = invalid_diagnostics),
    "'diagnostics' must be a dataframe obtained by diagnostics_dynamics()"
  )

  # 'divergences' must be logical (TRUE or FALSE)
  expect_error(
    plot_diagnostics(
      data = valid_data, diagnostics = valid_diagnostics,
      divergences = "yes"
    ),
    "'divergences' must be either 'TRUE' or 'FALSE'"
  )

  # 'max_treedepth' must be logical (TRUE or FALSE)
  expect_error(
    plot_diagnostics(
      data = valid_data, diagnostics = valid_diagnostics,
      max_treedepth = "yes"
    ),
    "'max_treedepth' must be either 'TRUE' or 'FALSE'"
  )

  # 'Rhat' must be logical (TRUE or FALSE)
  expect_error(
    plot_diagnostics(data = valid_data, diagnostics = valid_diagnostics, Rhat = 1),
    "'Rhat' must be either 'TRUE' or 'FALSE'"
  )

  # 'n_eff' must be logical (TRUE or FALSE)
  expect_error(
    plot_diagnostics(data = valid_data, diagnostics = valid_diagnostics, n_eff = NULL),
    "'n_eff' must be either 'TRUE' or 'FALSE'"
  )

  # Check if 'data' contains a column named 'time'
  invalid_data_no_time <- data.frame(metabolite = c("A", "B", "C"))
  expect_error(
    plot_diagnostics(data = invalid_data_no_time, diagnostics = valid_diagnostics),
    "'data' must contain a column named 'time'"
  )
})

test_that("plot_diagnostics:output_checks", {
  # generates plot when conditions are met
  plots <- plot_diagnostics(diagnostics = valid_diagnostics, data = valid_data, divergences = TRUE, max_treedepth = TRUE, Rhat = TRUE, n_eff = TRUE)

  # Check if each diagnostic condition is represented in the plots
  expect_true("divergences" %in% names(plots))
  expect_true("max_treedepth" %in% names(plots))
  expect_true("Rhat" %in% names(plots))
  expect_true("n_eff" %in% names(plots))

  # Test for no diagnostics (i.e., all FALSE)
  plots <- plot_diagnostics(
    diagnostics = valid_diagnostics, data = valid_data,
    divergences = FALSE, max_treedepth = FALSE,
    Rhat = FALSE, n_eff = FALSE
  )
  expect_equal(length(plots), 0) # No plots should be generated
})

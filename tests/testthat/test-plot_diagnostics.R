test_that("plot_diagnostics: input checks", {
  
  # Mock valid data for testing
  valid_data <- data.frame(time = c(1, 2, 3), metabolite = c("A", "B", "C"))
  valid_diagnostics <- data.frame(
    condition = c("dose1", "dose2"),
    treedepth_error = c(FALSE, TRUE),
    rhat.mu = c(1.01, 1.05),
    rhat.value = c(1.00, 1.02),
    neff.mu = c(500, 600),
    neff.value = c(1000, 1500)
  )
  
  # Invalid data mocks
  invalid_data <- list(time = c(1, 2, 3))
  invalid_diagnostics <- list(condition = c("dose1", "dose2"))
  
  # 'data' must be a dataframe or SummarizedExperiment object
  expect_error(plot_diagnostics(data = invalid_data, diagnostics = valid_diagnostics),
               "'data' must be a dataframe or colData of a SummarizedExperiment object")
  
  # 'diagnostics' must be a dataframe
  expect_error(plot_diagnostics(data = valid_data, diagnostics = invalid_diagnostics),
               "'diagnostics' must be a dataframe obtained by diagnostics_dynamics()")
  
  # 'divergences' must be logical (TRUE or FALSE)
  expect_error(plot_diagnostics(valid_data, valid_diagnostics, divergences = "yes"),
               "'divergences' must be either 'TRUE' or 'FALSE'")
  
  # 'max_treedepth' must be logical (TRUE or FALSE)
  expect_error(plot_diagnostics(valid_data, valid_diagnostics, max_treedepth = "yes"),
               "'max_treedepth' must be either 'TRUE' or 'FALSE'")
  
  # 'Rhat' must be logical (TRUE or FALSE)
  expect_error(plot_diagnostics(valid_data, valid_diagnostics, Rhat = 1),
               "'Rhat' must be either 'TRUE' or 'FALSE'")
  
  # 'n_eff' must be logical (TRUE or FALSE)
  expect_error(plot_diagnostics(valid_data, valid_diagnostics, n_eff = NULL),
               "'n_eff' must be either 'TRUE' or 'FALSE'")
  
  # Check if 'data' contains a column named 'time'
  invalid_data_no_time <- data.frame(metabolite = c("A", "B", "C"))
  expect_error(plot_diagnostics(invalid_data_no_time, valid_diagnostics),
               "'data' must contain a column named 'time'")
  
})
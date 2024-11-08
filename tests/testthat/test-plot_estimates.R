test_that("plot_estimates: input checks", {
  
  # Mock valid data for testing
  valid_data <- data.frame(time = c(1, 2, 3), metabolite = c("A", "B", "C"))
  valid_estimates <- data.frame(
    condition = c("dose1", "dose2"),
    metabolite.ID = c(1, 2),
    metabolite = c("A", "B"),
    higher = c(1.5, 2.5),
    lower = c(0.5, 1.5)
  )
  
  # Invalid data mocks
  invalid_data <- list(time = c(1, 2, 3))
  invalid_estimates <- list(condition = c("dose1", "dose2"))
  
  # 'data' must be a dataframe or SummarizedExperiment object
  expect_error(plot_estimates(estimates = valid_estimates, data = invalid_data),
               "'data' must be a dataframe or colData of a SummarizedExperiment object")
  
  # 'estimates' must be a dataframe
  expect_error(plot_estimates(invalid_estimates, valid_data),
               "'estimates' must be a dataframe obtained by estimates_dynamics()")
  
  # delta_t' must be logical 
  expect_error(plot_estimates(valid_estimates, valid_data, delta_t = "yes"),
               "'delta_t' must be either 'TRUE' or 'FALSE'")
  
  # 'dynamics' must be logical (TRUE or FALSE)
  expect_error(plot_estimates(valid_estimates, valid_data, dynamics = 1),
               "'dynamics' must be either 'TRUE' or 'FALSE'")
  
})

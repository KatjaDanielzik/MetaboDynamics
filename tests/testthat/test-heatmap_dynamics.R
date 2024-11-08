test_that("heatmap_dynamics: input checks", {
  
  # Mock data for testing
  clusters <- as.data.frame(cbind(
    cluster = 1:3,
    condition = c("A", "B", "C")
  ))
  datas <- as.data.frame(cbind(
    cluster = 1:3,
    condition = c("A", "B", "C")
  ))
  
  clusters_error <- list()
  data_error <- list()
  
  # Check 1: inputs are dataframes
  expect_error(heatmap_dynamics(estimates = data_error, clusters= clusters),
               "'estimates' must be a dataframe obtained by compare_dynamics()")
  
  expect_error(heatmap_dynamics(estimates = datas, clusters= clusters_error),
               "'clusters' must be a dataframe")

})
# Mock data for testing
clusters <- list(A=list(data=as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))))
data <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))

test_that("heatmap_dynamics: input checks", {
  clusters_error <- data.frame(
  )
  data_error <- list()

  # Check 1: inputs are dataframes
  expect_error(
    heatmap_dynamics(estimates = data_error, data = clusters),
    "'estimates' must be a dataframe obtained by compare_dynamics()"
  )

  expect_error(
    heatmap_dynamics(estimates = data, data = clusters_error),
    "'data' must be list or a SummarizedExperiment object"
  )
})

test_that("heatmap_dynamics:output_checks", {
  plot <- heatmap_dynamics(estimates = data, data = clusters )
  plot_title <- plot$labels$title
  expect_equal(plot_title, "similarity of dynamics in clusters")

  plot_data <- plot$data
  expect_true(all(plot_data$size > 0)) # size should be positive (1 / uncertainty)
  expect_true(all(plot_data$col > 0)) # color should be positive (1 / mu_mean)
})

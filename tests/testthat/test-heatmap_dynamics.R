# Mock data for testing
clusters <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))
datas <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))

test_that("heatmap_dynamics: input checks", {
  clusters_error <- list()
  data_error <- list()

  # Check 1: inputs are dataframes
  expect_error(
    heatmap_dynamics(estimates = data_error, data = clusters),
    "'estimates' must be a dataframe obtained by compare_dynamics()"
  )

  expect_error(
    heatmap_dynamics(estimates = datas, data = clusters_error),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )
})

test_that("heatmap_dynamics:output_checks", {
  plot <- heatmap_dynamics(estimates = datas, data = clusters)
  plot_title <- plot$labels$title
  expect_equal(plot_title, "similarity of dynamics in clusters")

  # check the behavior with missing 'mu_mean' or CrI values
  # Create mock data with NA values
  estimates_with_na <- datas
  estimates_with_na$mu_mean[1] <- NA
  estimates_with_na$`97.5%`[1] <- NA
  estimates_with_na$`2.5%`[1] <- NA
  plot <- heatmap_dynamics(estimates = estimates_with_na, data = clusters)
  # Check if the plot can still be generated (no errors)
  expect_s3_class(plot, "gg") # Ensure the output is still a ggplot object

  # Check that 'mu_mean' and uncertainty (CrI) are correctly reflected in plot size and color
  plot <- heatmap_dynamics(estimates = datas, data = clusters)
  # Check that the size and color of points are being set as expected
  plot_data <- plot$data
  expect_true(all(plot_data$size > 0)) # size should be positive (1 / uncertainty)
  expect_true(all(plot_data$col > 0)) # color should be positive (1 / mu_mean)


  # Check for empty input
  empty_estimates <- data.frame(
    parameter = character(0), cluster_a = integer(0),
    cluster_b = integer(0), mu_mean = numeric(0),
    `97.5%` = numeric(0), `2.5%` = numeric(0)
  )
  empty_clusters <- data.frame(cluster = integer(0), condition = character(0))
  plot <- heatmap_dynamics(estimates = empty_estimates, data = empty_clusters)
  # Check that the plot is still returned, even if empty data is provided
  expect_s3_class(plot, "gg")
})

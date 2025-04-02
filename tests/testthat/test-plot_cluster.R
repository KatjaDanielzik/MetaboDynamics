test_that("plot_cluster: input checks", {
  # Function returns an error when data is not a list of dataframes or SummarizedExperiment object
  test_data <- data.frame(a = c(1, 2, 3))
  expect_error(
    plot_cluster(test_data),
    "'data' must be a list or a SummarizedExperiment object
         obtained by function cluster_dynamics()"
  )
})

test_that("plot_cluster: output checks", {
  # Function returns a list with the correct number of elements when data is a list of dataframes obtained by function cluster_dynamics()

  test_data <- list(list(
    dynamics = c("X1", "X2"),
    data = data.frame(
      metabolite = c("ATP", "L-Alanine", "GDP"),
      "1" = c(1, 2, 3),
      "2" = c(1, 2, 3),
      cluster = c(1, 2, 1),
      condition = "A"
    ),
    tree = hclust(dist(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)))
  ))
  result <- plot_cluster(test_data)
  expect_equal(length(result), 2)
  expect_true(inherits(result[[1]]$dendrogram, "recordedplot"))
  expect_true(inherits(result[[1]]$PCA_plot, "ggplot"))
  expect_true(inherits(result$lineplot, "ggplot"))
})

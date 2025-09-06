data("longitudinalMetabolomics")
dummy_data <- longitudinalMetabolomics 
dummy_cluster <- metadata(dummy_data)[["cluster"]]

# Define tests
test_that("plot_cluster:input_checks", {
  # Test that the function throws an error if data is not a list or a SummarizedExperiment object
  expect_error(
    plot_cluster(data = "not a list or a SummarizedExperiment object"),
    "'data' must be a list or a SummarizedExperiment object obtained by function cluster_dynamics()"
  )
  
test_that("plot_cluster:output_checks", {
  # Test that the function returns a list with the expected plot structure
  results <- plot_cluster(data = dummy_cluster)
  expect_type(results, "list")
  expected_elements <- c("trees", "clusterplots", "lineplots", "patchwork", "cluster_order","cluster_heights")
  expect_true(all(expected_elements %in% names(results)))
  
  # Test specific components of the output
  expect_true(all(sapply(results$trees, function(x) inherits(x, "gg"))))
  expect_true(all(sapply(results$clusterplots, function(x) inherits(x, "gg"))))
  expect_true(all(sapply(results$lineplots, function(x) inherits(x, "patchwork"))))
  expect_true(all(sapply(results$patchwork, function(x) inherits(x, "patchwork"))))
  
  # Test that the function processes SummarizedExperiment objects correctly
  results_se <- plot_cluster(data = dummy_data)
  expect_type(results_se, "list")
  expect_true(all(expected_elements %in% names(results_se)))
})
})
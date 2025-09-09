# Define a dummy data frame
dummy_df <- data.frame(
  condition = rep(c("CondA", "CondB"), each = 4),
  cluster = rep(1:2, times = 4),
  middle_hierarchy = rep(c("Mod1", "Mod2"), times = 4),
  OvE_gen = runif(8, 0.5, 1.5),
  OvE_gen_median = runif(8, 0.5, 1.5),
  OvE_gen_lower = runif(8, 0.1, 0.8),
  OvE_gen_higher = runif(8, 1.2, 2)
)

# Define a dummy SummarizedExperiment object
dummy_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = matrix(rnorm(20), ncol = 4)),
  colData = data.frame(x = 1:4),
  metadata = list(
    ORA_middle_hierarchy = dummy_df
  )
)

# Define a dummy result from plot_cluster() for patchwork testing
dummy_cluster_plot_result <- list(
  cluster_order = list(
    CondA = c("1", "2"),
    CondB = c("1", "2")
  ),
  cluster_heights = list(
    CondA = c(2, 3),
    CondB = c(4, 1)
  )
)

# Define tests
test_that("plot_ORA:input_checks", {
  # Test that the function throws an error if data is not a data.frame or SummarizedExperiment
  expect_error(
    plot_ORA(data = "not a data.frame or SummarizedExperiment"),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )

  # Test that the function throws an error if patchwork is not logical
  expect_error(
    plot_ORA(data = dummy_df, patchwork = "not logical"),
    "'patchwork' must be either 'TRUE' or 'FALSE'"
  )

  # Test that the function throws an error if patchwork is TRUE but plot_cluster is missing or not a list
  expect_error(
    plot_ORA(data = dummy_df, patchwork = TRUE),
    "if 'patchwork is TRUE, plot_cluster must be provided"
  )
  expect_error(
    plot_ORA(data = dummy_df, patchwork = TRUE, plot_cluster = "not a list"),
    "if 'patchwork is TRUE, plot_cluster must be provided"
  )
})

test_that("plot_ORA:output_checks", {
  # Test that the function returns a list with the expected elements
  results <- plot_ORA(data = dummy_df, patchwork = FALSE)
  expect_type(results, "list")
  expect_named(results, c("plot", "ora_patchwork"))
  expect_s3_class(results$plot, "gg")
  expect_length(results$ora_patchwork, 0)

  # Test that the function returns a list with the expected elements when patchwork is TRUE
  results_patch <- plot_ORA(data = dummy_df, patchwork = TRUE, plot_cluster = dummy_cluster_plot_result)
  expect_type(results_patch, "list")
  expect_named(results_patch, c("plot", "ora_patchwork"))
  expect_s3_class(results_patch$plot, "gg")
  expect_true(is.list(results_patch$ora_patchwork))
})

test_that("plot_ORA:output_checks_with_se", {
  # Test that the function returns a list with the expected elements when using a SummarizedExperiment object
  results_se <- plot_ORA(data = dummy_se, patchwork = FALSE)
  expect_type(results_se, "list")
  expect_named(results_se, c("plot", "ora_patchwork"))
  expect_s3_class(results_se$plot, "gg")
  expect_length(results_se$ora_patchwork, 0)
})

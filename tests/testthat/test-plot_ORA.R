test_that("plot_ORA: input checks", {
  invalid_ORA <- list(
    OvE_gen = c(1.5, 0.8),
    module_name = c("module1", "module2")
  )

  # ORA is a dataframe
  expect_error(
    plot_ORA(invalid_ORA),
    "'data' must be a dataframe or a SummarizedExperiment object"
  )
})

test_that("plot_ORA:output_checks", {
  # OvE have to be positive
  ORA_example <- data.frame(
    cluster = rep(c("Cluster1", "Cluster2"), each = 5),
    condition = rep(c("Condition1", "Condition2"), each = 5),
    module_name = paste("Module", 1:10),
    OvE_gen = abs(rnorm(10)),
    OvE_gen_lower = abs(abs(rnorm(10)) - 1),
    OvE_gen_higher = abs(rnorm(10)) + 1,
    OvE_gen_median = abs(rnorm(10))
  )

  plot <- plot_ORA(ORA_example)
  expect_s3_class(plot, "gg")
})

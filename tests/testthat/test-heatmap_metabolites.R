# Mock data for testing
clusters <- list(A=list(data=as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))))

data <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))

test_that("heatmap_metabolites: input checks", {
  clusters_error <- data.frame()
  data_error <- list()

  # Check 1: inputs are dataframes
  expect_error(
    heatmap_metabolites(distances = data_error, data = clusters),
    "'distances' must be a dataframe obtained by compare_metabolites()"
  )

  expect_error(
    heatmap_metabolites(distances = data, data = clusters_error)
  )
})


test_that("heatmap_metabolites:output_checks", {
  # returns ggplot objext
  plot <- heatmap_metabolites(distances = data, data = clusters)
  expect_s3_class(plot, "gg")

  # empty input without crashing
  empty_distances <- data.frame(Var1 = integer(0), Var2 = integer(0), Jaccard = numeric(0))
  empty_clusters <- data.frame(cluster = integer(0), condition = character(0))
})

# Mock data for testing
clusters <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))
datas <- as.data.frame(cbind(
  cluster = 1:3,
  condition = c("A", "B", "C")
))

test_that("heatmap_metabolites: input checks", {
  clusters_error <- list()
  data_error <- list()

  # Check 1: inputs are dataframes
  expect_error(
    heatmap_metabolites(distances = data_error, data = clusters),
    "'distances' must be a dataframe obtained by compare_metabolites()"
  )

  expect_error(
    heatmap_metabolites(distances = datas, data = clusters_error),
    "'data' must be a dataframe"
  )
})


test_that("heatmap_metabolites:output_checks", {
  # returns ggplot objext
  plot <- heatmap_metabolites(distances = datas, data = clusters)
  expect_s3_class(plot, "gg")

  # empty input without crashing
  empty_distances <- data.frame(Var1 = integer(0), Var2 = integer(0), Jaccard = numeric(0))
  empty_clusters <- data.frame(cluster = integer(0), condition = character(0))

  plot <- heatmap_metabolites(distances = empty_distances, data = empty_clusters)
  # Check that the plot is still returned, even if empty data is provided
  expect_s3_class(plot, "gg")

  distances_with_na <- datas
  distances_with_na$Jaccard[1] <- NA

  # NAs in Jaccard column
  plot <- heatmap_metabolites(distances = distances_with_na, data = clusters)
  # Check that the plot is still returned even if there are NA values
  expect_s3_class(plot, "gg")
})

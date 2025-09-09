test_that("compare_dynamics: input checks", {
  error_clusters <- list(A = data.frame(
    condi = c(rep("A", 5), rep("B", 5)),
    clust = c(rep("C1", 5), rep("C2", 5)),
    mu1_mean = rnorm(10), mu2_mean = rnorm(10)
  ))

  expect_error(
    compare_dynamics(
      data = NULL
    )
  )
  expect_error(
    compare_dynamics(
      data = error_clusters
    )
  )
})

test_that("compare_dynamics:output_checks", {
  # Prepare valid dummy data
  dummy_clusters <- list(A = list(data = data.frame(
    metabolite = rep("A",10),
    condition = rep(c("A", "B"), each = 5),
    cluster = rep(c("1", "2"), times = 5),
    mu1_mean = rnorm(10), mu2_mean = rnorm(10)
  )))

  result <- compare_dynamics(data = dummy_clusters, cores = 1)

  #  Output is a list with expected names
  expect_type(result, "list")
  expect_named(result, c("distances", "fit", "estimates"))

  # 'distances' is a list of numeric vectors
  expect_type(result$distances, "list")
  expect_true(all(sapply(result$distances, is.numeric)))

  # 'fit' is a stanfit object
  expect_s4_class(result[["fit"]], "stanfit")

  # 'estimates' is a data frame with expected columns
  expect_s3_class(result$estimates, "data.frame")
  expect_true(all(c("comparison", "mu_mean", "sigma_mean", "cluster_a", "cluster_b") %in% colnames(result$estimates)))

  # Number of comparisons
  n_expected <- choose(4, 2) # For 4 unique combinations (2 conditions à 2 clusters), 6 pairwise comparisons
  expect_equal(length(result$distances), n_expected)
})

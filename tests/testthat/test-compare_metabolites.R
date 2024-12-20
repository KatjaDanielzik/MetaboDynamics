test_that("compare_metabolites: input checks", {
  mock_clusters <- data.frame(
    metabolite = c(rep("A", 5), rep("B", 5)),
    cluster = c(rep("C1", 5), rep("C2", 5))
  )

  error_clusters <- data.frame(
    met = c(rep("A", 5), rep("B", 5)),
    clust = c(rep("C1", 5), rep("C2", 5))
  )

  expect_error(
    compare_metabolites(data = NULL, metabolite = "metabolite"),
    "'data' must be a dataframe"
  )
  expect_error(
    compare_metabolites(data = mock_clusters, metabolite = NULL),
    "'metabolite' must be a character vector"
  )
  expect_error(
    compare_metabolites(
      data = error_clusters,
      metabolite = "metabolite"
    ),
    "'data' must contain a column named 'cluster'"
  )
  expect_error(
    compare_metabolites(data = mock_clusters, metabolite = "met"),
    "'data' must contain a column containing
         metabolite names as specified with metabolite= "
  )
  expect_error(compare_metabolites(clusters = error_clusters, metabolite = "met"))
})

test_that("compare_metabolites:output_checks", {
  # Prepare valid dummy data
  dummy_clusters <- data.frame(
    condition = rep(c("A", "B"), each = 3),
    cluster = rep(c("1", "2"), times = 3),
    metabolite = c("met1", "met2", "met3", "met1", "met4", "met5")
  )

  result <- compare_metabolites(data = dummy_clusters, metabolite = "metabolite")

  # Output is a data frame
  expect_s3_class(result, "data.frame")

  # Output contains expected columns
  expected_columns <- c("cluster_a", "cluster_b", "Jaccard")
  expect_true(all(expected_columns %in% colnames(result)))

  # Jaccard column values are numeric
  expect_type(result$Jaccard, "double")

  # Jaccard index is between 0 and 1
  expect_true(all(result$Jaccard >= 0 & result$Jaccard <= 1))

  # Number of comparisons matches expected (nrow(x) * nrow(x))
  unique_clusters <- unique(dummy_clusters[, c("condition", "cluster")])
  # Generate expected comparisons (all pairs between group_a and group_b)
  n_combinations <- nrow(unique_clusters)
  # Generate upper triangle of pairwise combinations (i.e., no self-comparisons)
  comparison_pairs <- combn(seq_len(n_combinations), 2)
  expect_equal(nrow(result), ncol(comparison_pairs))
})

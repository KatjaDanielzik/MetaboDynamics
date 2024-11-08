test_that("compare_dynamics: input checks", {
  mock_clusters <- data.frame(
    condition = c(rep("A", 5), rep("B", 5)),
    cluster = c(rep("C1", 5), rep("C2", 5)),
    mu1_mean = rnorm(10), mu2_mean = rnorm(10)
  )
  
  error_clusters <-  data.frame(
    condi = c(rep("A", 5), rep("B", 5)),
    clust = c(rep("C1", 5), rep("C2", 5)),
    mu1_mean = rnorm(10), mu2_mean = rnorm(10)
  )
  
  expect_error(compare_dynamics(clusters = NULL, 
                                dynamics = c("mu1_mean","mu2_mean","mu3_mean","mu4_mean")), 
               "'clusters' must be a dataframe")
  expect_error(compare_dynamics(clusters = mock_clusters, dynamics = NULL), 
               "'dynamics' must be a character vector")
  expect_error(compare_dynamics(clusters = error_clusters, 
                                dynamics = c("mu1_mean","mu2_mean","mu3_mean","mu4_mean")), 
  "'clusters' must contain 'condition' and 'cluster' columns")
  expect_error(compare_dynamics(clusters = mock_clusters, dynamics = c("not_present")), 
               "All specified 'dynamics' columns must exist in `clusters` dataframe")
})


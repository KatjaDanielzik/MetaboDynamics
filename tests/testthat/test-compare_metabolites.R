test_that("compare_dynamics: input checks", {
  mock_clusters <- data.frame(
    metabolite = c(rep("A", 5), rep("B", 5)),
    cluster = c(rep("C1", 5), rep("C2", 5))
  )
  
  error_clusters <-  data.frame(
    met = c(rep("A", 5), rep("B", 5)),
    clust = c(rep("C1", 5), rep("C2", 5))
  )
  
  expect_error(compare_metabolites(clusters = NULL, metabolite="metabolite"), 
               "'clusters' must be a dataframe")
  expect_error(compare_metabolites(clusters = mock_clusters, metabolite = NULL), 
               "'metabolite' must be a character vector")
  expect_error(compare_metabolites(clusters = error_clusters, 
                                metabolite = "metabolite"), 
               "'clusters' must contain a column named 'cluster'")
  expect_error(compare_metabolites(clusters = mock_clusters, metabolite="met"), 
               "'clusters' must contain a column containing
         metabolite names as specified with metabolite= ")
  expect_error(compare_metabolites(clusters = error_clusters, metabolite="met"))
})

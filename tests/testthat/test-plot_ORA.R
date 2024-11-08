test_that("plot_ORA: input checks", {
  
 invalid_ORA <- list(
    OvE_gen = c(1.5, 0.8),
    module_name = c("module1", "module2")
  )
  
  # ORA is a dataframe
  expect_error(
    plot_ORA(invalid_ORA),
    "'ORA' must be a dataframe obtained by ORA_hypergeometric()")
})
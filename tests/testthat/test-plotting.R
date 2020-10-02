context("Test Plotting Functions")
library("BAGEL")

test_that(desc = "Testing Exposures Plotting", {
  result <- readRDS(system.file("testdata", "res.rds", package = "BAGEL"))
  p <- plot_exposures(result)
  expect_true(ggplot2::is.ggplot(p))
})

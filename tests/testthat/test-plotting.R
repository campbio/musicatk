context("Test Plotting Functions")
library("musicatk")

test_that(desc = "Testing Exposures Plotting", {
  result <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
  p <- plot_exposure_bar(result)
  expect_true(ggplot2::is.ggplot(p))
})

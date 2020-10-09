context("Test Plotting Functions")
library("musicatk")

test_that(desc = "Testing Exposures Plotting", {
  result <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
  p <- plot_exposures(result, plot_type = "bar")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_exposures(result, plot_type = "violin")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_exposures(result, plot_type = "box")
  expect_true(ggplot2::is.ggplot(p))
})

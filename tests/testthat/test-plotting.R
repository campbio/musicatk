context("Test Plotting Functions")
library("musicatk")

test_that(desc = "Testing Exposures Plotting", {
  data(res)
  p <- plot_exposures(res, model_name = "res", plot_type = "bar")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_exposures(res, model_name = "res", plot_type = "violin")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_exposures(res, model_name = "res", plot_type = "box")
  expect_true(ggplot2::is.ggplot(p))
})

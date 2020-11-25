context("Test Heatmap plotting function")
library("musicatk")

test_that(desc = "Testing Exposures Heatmap Plotting", {
  result <- readRDS(system.file("testdata", "res_annot.rds", package = "musicatk"))
  p <- plot_heatmap(result)
  expect_true(!is.null(p))
  p <- plot_heatmap(result, proportional = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(result, proportional = TRUE,col_annot = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(result, col_annot = TRUE )
  expect_true(!is.null(p))
  #p <- plot_heatmap(result, col_samps)
  #expect_true(!is.null(p))
  
  
})
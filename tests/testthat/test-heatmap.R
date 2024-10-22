context("Test Heatmap plotting function")
library("musicatk")

test_that(desc = "Testing Exposures Heatmap Plotting", {
  data(res_annot)
  p <- plot_heatmap(res_annot, "res_annot")
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", proportional = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", proportional = TRUE,scale = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_column_names = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_row_names = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_column_names = TRUE, show_row_names = TRUE)
  expect_true(!is.null(p))
  
})

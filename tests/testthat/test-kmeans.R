context("Test Kmeans and plotting functions")
library(musicatk)

test_that(desc = "Test plotting function for selecting k", {
  data(res_annot)
  p <- k_select(res_annot, method = "elbow", n = 6)
  expect_true(ggplot2::is.ggplot(p))
  p <- k_select(res_annot, method = "silhouette", n = 6)
  expect_true(ggplot2::is.ggplot(p))
  p <- k_select(res_annot, method = "gap", n =6)
  expect_true(ggplot2::is.ggplot(p))
})

test_that(desc = "Test kmeans visualization function", {
  data(res_annot)
  k_out <- cluster_kmeans(result = res_annot, centers = 2)
  expect_equal(nrow(k_out), 7)
  expect_error(plot_kmeans(result = res_annot, clusters = k_out, group = "signature"), "UMAP not found")
  expect_error(plot_kmeans(result = res_annot, clusters = k_out, group = "annotation"), "UMAP not found")
  expect_error(plot_kmeans(result = res_annot, clusters = k_out, group = "none"), "UMAP not found")
  create_umap(res_annot)
  test_res <- res_annot
  test_res@musica@sample_annotations <- test_res@musica@sample_annotations[,-2]
  expect_error(plot_kmeans(result = test_res, clusters = k_out, group = "annotation"), "Sample annotation not found")
  p <- plot_kmeans(result = res_annot, clusters = k_out, group = "signature")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_kmeans(result = res_annot, clusters = k_out, group = "annotation")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_kmeans(result = res_annot, clusters = k_out, group = "none")
  expect_true(ggplot2::is.ggplot(p))
})
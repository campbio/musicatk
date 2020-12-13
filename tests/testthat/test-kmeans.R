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
  clust_out <- cluster_exposure(result = res_annot, nclust = 2)
  expect_equal(nrow(clust_out), 7)
  expect_error(plot_kmeans(result = res_annot, clusters = clust_out, group = "signature"), "UMAP not found")
  expect_error(plot_kmeans(result = res_annot, clusters = clust_out, group = "annotation"), "UMAP not found")
  expect_error(plot_kmeans(result = res_annot, clusters = clust_out, group = "none"), "UMAP not found")
  create_umap(res_annot)
  test_res <- res_annot
  test_res@musica@sample_annotations <- test_res@musica@sample_annotations[,-2]
  expect_error(plot_kmeans(result = test_res, clusters = clust_out, group = "annotation"), "Sample annotation not found")
  expect_error(plot_kmeans(result = res_annot, clusters = clust_out, group = "annotation", annotation = "cancer"), "invalid annotation column name")
  p <- plot_kmeans(result = res_annot, clusters = clust_out, group = "signature")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_kmeans(result = res_annot, clusters = clust_out, group = "annotation", annotation = "Tumor_Subtypes")
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_kmeans(result = res_annot, clusters = clust_out, group = "none")
  expect_true(ggplot2::is.ggplot(p))
})
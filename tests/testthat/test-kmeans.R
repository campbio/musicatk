context("Test Kmeans and plotting functions")
library(musicatk)

test_that(desc = "Test plotting function for selecting k", {
  data(res_annot)
  p <- k_select(res_annot, "res_annot", method = "wss", n = 6)
  expect_true(ggplot2::is.ggplot(p))
  p <- k_select(res_annot, "res_annot", method = "silhouette", n = 6)
  expect_true(ggplot2::is.ggplot(p))
  p <- k_select(res_annot, "res_annot", method = "gap_stat", n = 6)
  expect_true(ggplot2::is.ggplot(p))
})

test_that(desc = "Test kmeans visualization function", {
  data(res_annot)
  expect_error(cluster_exposure(musica = res_annot, model_name = "res_annot",
                                nclust = 2, method = "hkmeans",
                                dis.method = "cosine"),
               "For hkmeans clustering")
  expect_error(cluster_exposure(musica = res_annot, model_name = "res_annot",
                                nclust = 2, method = "clara",
                                dis.method = "cosine"), "For clara clustering")
  clust_out <- cluster_exposure(musica = res_annot, model_name = "res_annot",
                                nclust = 2, iter.max = 15)
  expect_equal(nrow(clust_out), 7)
  expect_error(plot_cluster(musica = res_annot, model_name = "res_annot",
                            clusters = clust_out, group = "signature"),
               "UMAP not found")
  expect_error(plot_cluster(musica = res_annot, model_name = "res_annot",
                            clusters = clust_out, group = "annotation"),
               "UMAP not found")
  expect_error(plot_cluster(musica = res_annot, model_name = "res_annot",
                            clusters = clust_out, group = "none"),
               "UMAP not found")
  create_umap(res_annot, model_name = "res_annot")
  test_res <- res_annot
  test_res@sample_annotations <- as.data.frame(samp_annot(res_annot)[, -2])
  colnames(test_res@sample_annotations) <- "Sample"
  expect_error(plot_cluster(musica = test_res, model_name = "res_annot",
                            clusters = clust_out, group = "annotation"),
               "Sample annotation not found")
  expect_error(plot_cluster(musica = res_annot, model_name = "res_annot",
                            clusters = clust_out, group = "annotation",
                            annotation = "cancer"),
               "invalid annotation column name")
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "signature", plotly = FALSE)
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "annotation",
                    annotation = "Tumor_Subtypes", plotly = FALSE)
  expect_true(ggplot2::is.ggplot(p))
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "none", plotly = FALSE)
  expect_true(ggplot2::is.ggplot(p))
})

context("Test Kmeans and plotting functions")
library(musicatk)

test_that(desc = "Test plotting function for selecting k", {
  data(res_annot)
  p <- k_select(res_annot, "res_annot", method = "wss", n = 6)
  expect_s3_class(p, "ggplot")
  p <- k_select(res_annot, "res_annot", method = "silhouette", n = 6)
  expect_s3_class(p, "ggplot")
  p <- k_select(res_annot, "res_annot", method = "gap_stat", n = 6)
  expect_s3_class(p, "ggplot")
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
  expect_equal(nrow(clust_out), 12)
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
  expect_s3_class(p, "ggplot")
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "annotation",
                    annotation = "Tumor_Subtypes", plotly = FALSE)
  expect_s3_class(p, "ggplot")
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "none", plotly = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that(desc = "plot_cluster uses the chosen annotation, matched by name", {
  data(res_annot)
  samp_annot(res_annot, "test") <- rep(c("A", "B", "C"), 4)
  clust_out <- cluster_exposure(musica = res_annot, model_name = "res_annot",
                                nclust = 2, iter.max = 15)
  create_umap(res_annot, model_name = "res_annot")
  p <- plot_cluster(musica = res_annot, model_name = "res_annot",
                    clusters = clust_out, group = "annotation",
                    annotation = "test", plotly = FALSE)
  umap_samples <- rownames(umap(res_annot, "result", "SBS96", "res_annot"))
  annot <- samp_annot(res_annot)
  expected <- annot$test[match(umap_samples, as.character(annot$Samples))]
  expect_equal(as.character(p$data$annotation), expected)
})

# Proportional exposures of res_annot, samples in rows, as clustering uses
.clust_fixture <- function() {
  data(res_annot, envir = environment())
  e <- exposures(res_annot, "result", "SBS96", "res_annot")
  list(musica = res_annot, expos = t(sweep(e, 2, colSums(e), FUN = "/")))
}

.clust_diss <- function(expos, method) {
  d <- philentropy::distance(expos,
    method = method, use.row.names = TRUE,
    as.dist.obj = TRUE, mute.message = TRUE
  )
  if (method == "cosine") {
    d <- 1 - d
  }
  d
}

test_that(desc = "hclust uses the chosen dissimilarity (#173)", {
  fx <- .clust_fixture()
  # k where devel's distances-of-distances gave different clusters
  cases <- list(euclidean = 6, manhattan = 6, cosine = 2)
  for (m in names(cases)) {
    expected <- stats::cutree(
      stats::hclust(.clust_diss(fx$expos, m), method = "ward.D"),
      k = cases[[m]]
    )
    observed <- cluster_exposure(fx$musica, "res_annot",
      nclust = cases[[m]], method = "hclust", dis.method = m
    )
    expect_equal(rownames(observed), names(expected))
    expect_equal(observed$cluster, unname(expected))
  }
})

test_that(desc = "cosine dis.method is 1 - cosine similarity", {
  fx <- .clust_fixture()
  expected <- cluster::pam(.clust_diss(fx$expos, "cosine"),
    k = 2, diss = TRUE, cluster.only = TRUE
  )
  observed <- cluster_exposure(fx$musica, "res_annot",
    nclust = 2, method = "pam", dis.method = "cosine"
  )
  expect_equal(observed$cluster, unname(expected))
})

test_that(desc = "kmeans runs k-means; fanny runs fuzzy clustering (#34)", {
  fx <- .clust_fixture()
  expect_warning(
    observed <- cluster_exposure(fx$musica, "res_annot",
      nclust = 4, method = "kmeans"
    ),
    NA
  )
  expected <- withr::with_seed(1, stats::kmeans(fx$expos,
    centers = 4,
    iter.max = 10, nstart = 10
  ))
  expect_equal(rownames(observed), names(expected$cluster))
  expect_equal(observed$cluster, unname(expected$cluster))
  # Seeded: the same clusters each call, and the caller's RNG is untouched
  set.seed(42)
  before <- .Random.seed
  again <- cluster_exposure(fx$musica, "res_annot",
    nclust = 4, method = "kmeans"
  )
  expect_true(identical(before, .Random.seed))
  expect_equal(again, observed)

  expect_error(
    cluster_exposure(fx$musica, "res_annot",
      nclust = 2, method = "kmeans", dis.method = "cosine"
    ),
    "euclidean"
  )

  expected <- cluster::fanny(.clust_diss(fx$expos, "euclidean"),
    k = 2, diss = TRUE, maxit = 500, tol = 1e-15
  )
  expect_warning(
    observed <- cluster_exposure(fx$musica, "res_annot",
      nclust = 2, method = "fanny"
    ),
    NA
  )
  expect_equal(observed$cluster, unname(expected$clustering))
})

test_that(desc = "k_select lowers n to the number of samples - 1 (#185)", {
  fx <- .clust_fixture()
  .keep_samples <- function(musica, n) {
    e <- exposures(musica, "result", "SBS96", "res_annot")
    exposures(musica, "result", "SBS96", "res_annot") <-
      e[, seq_len(n), drop = FALSE]
    musica
  }
  ten <- .keep_samples(fx$musica, 10)
  for (clust in c("kmeans", "pam")) {
    expect_message(
      p <- k_select(ten, "res_annot", method = "wss", clust.method = clust),
      "9"
    )
    expect_s3_class(p, "ggplot")
  }
  expect_error(
    k_select(.keep_samples(fx$musica, 2), "res_annot"),
    "at least three samples"
  )
  expect_error(k_select(fx$musica, "res_annot", n = NA), "'n' must be")
})

test_that(desc = "plot_cluster gives no ggplot2 deprecation warnings", {
  data(res_annot)
  withr::local_options(lifecycle_verbosity = "warning")
  suppressMessages(create_umap(res_annot, "res_annot"))
  clust_out <- cluster_exposure(res_annot, "res_annot", nclust = 2)
  # expect_no_warning() lets lifecycle deprecation warnings through
  warnings <- character()
  withCallingHandlers({
    plot_cluster(res_annot, "res_annot", clusters = clust_out,
                 group = "none")
    plot_cluster(res_annot, "res_annot", clusters = clust_out,
                 group = "annotation", annotation = "Tumor_Subtypes")
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("deprecated", warnings)))
})

test_that(desc = "cluster_exposure names samples with zero exposure", {
  data(res_annot)
  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  exposure_matrix[, 4] <- 0
  exposures(res_annot, "result", "SBS96", "res_annot") <- exposure_matrix
  expect_error(cluster_exposure(res_annot, "res_annot", nclust = 2),
               regexp = colnames(exposure_matrix)[4])
})

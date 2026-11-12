context("Test UMAP functions")
library("musicatk")

test_that(desc = "plot_umap before create_umap says to run create_umap", {
  data(res_annot)
  expect_error(plot_umap(res_annot, "res_annot"), regexp = "create_umap")
})

test_that(desc = "create_umap names samples with zero exposure", {
  data(res_annot)
  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  exposure_matrix[, 3] <- 0
  exposures(res_annot, "result", "SBS96", "res_annot") <- exposure_matrix
  expect_error(create_umap(res_annot, "res_annot"),
               regexp = colnames(exposure_matrix)[3])
})

# Colours of the points in each panel of a signature UMAP
.umap_point_colours <- function(p) {
  ggplot2::layer_data(p, 1)$colour
}

test_that(desc = "plot_umap colour limits include the top samples", {
  data(res_annot)
  suppressMessages(create_umap(res_annot, "res_annot"))
  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  # Each signature's maximum, 4549, rounds down to 4500 with signif()
  exposure_matrix <- round(exposure_matrix /
                             apply(exposure_matrix, 1, max) * 4549)
  exposures(res_annot, "result", "SBS96", "res_annot") <- exposure_matrix

  p <- plot_umap(res_annot, "res_annot", proportional = FALSE)
  expect_false(any(.umap_point_colours(p) == "grey50"))

  panels <- NULL
  local_mocked_bindings(
    grid.arrange = function(...) {
      panels <<- list(...)
      invisible(NULL)
    },
    .package = "gridExtra"
  )
  plot_umap(res_annot, "res_annot", proportional = FALSE,
            same_scale = FALSE)
  colours <- unlist(lapply(panels, .umap_point_colours))
  expect_false(any(colours == "grey50"))
})

test_that(desc = "plot_umap centroid layer maps y to UMAP_2", {
  data(res_annot)
  suppressMessages(create_umap(res_annot, "res_annot"))
  p <- plot_umap(res_annot, "res_annot", color_by = "annotation",
                 annotation = "Tumor_Subtypes", add_annotation_labels = TRUE)
  centroid <- p$layers[[2]]
  expect_identical(rlang::as_label(centroid$mapping$y), "UMAP_2")
})

test_that(desc = "plot_umap(color_by = 'cluster') needs clust, not annotation", {
  data(res_annot)
  suppressMessages(create_umap(res_annot, "res_annot"))
  clust <- cluster_exposure(res_annot, "res_annot", nclust = 2)
  p <- plot_umap(res_annot, "res_annot", color_by = "cluster", clust = clust)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$annotation,
               as.character(clust[p$data$sample, "cluster"]))
  expect_error(plot_umap(res_annot, "res_annot", color_by = "cluster"),
               regexp = "cluster_exposure")
})

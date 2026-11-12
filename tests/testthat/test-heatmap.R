context("Test Heatmap plotting function")
library("musicatk")

test_that(desc = "Testing Exposures Heatmap Plotting", {
  data(res_annot)
  p <- plot_heatmap(res_annot, "res_annot")
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", proportional = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", proportional = TRUE, scale = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_column_names = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_row_names = TRUE)
  expect_true(!is.null(p))
  p <- plot_heatmap(res_annot, "res_annot", show_column_names = TRUE, 
                    show_row_names = TRUE)
  expect_true(!is.null(p))
})

test_that(desc = "Heatmap annotations are matched to exposures by name", {
  data(res_annot)
  samp_annot(res_annot, "test") <- rep(c("A", "B", "C"), 4)
  # A model holding a subset of the samples
  results <- result_list(res_annot)
  model <- results[["result"]]@modality[["SBS96"]][["res_annot"]]
  model@exposures <- model@exposures[, -2]
  results[["result"]]@modality[["SBS96"]][["res_annot"]] <- model
  result_list(res_annot) <- results

  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  annot <- samp_annot(res_annot)
  expected <- annot$test[
    match(colnames(exposure_matrix), as.character(annot$Samples))
  ]

  p <- plot_heatmap(res_annot, "res_annot", annotation = "test")
  observed <- p@top_annotation@anno_list$test@fun@var_env$value
  expect_equal(unname(observed), expected)
})

test_that(desc = "Heatmap keeps annotation names that are not syntactic", {
  data(res_annot)
  subtypes <- samp_annot(res_annot)$Tumor_Subtypes
  samp_annot(res_annot, "Tumor subtype") <- subtypes
  p <- plot_heatmap(res_annot, "res_annot", annotation = "Tumor subtype")
  expect_identical(names(p@top_annotation@anno_list), "Tumor subtype")
})

test_that(desc = "subset_tumor keeps samples with an exact annotation match", {
  data(res_annot)
  annot <- samp_annot(res_annot)
  lung <- as.character(annot$Samples[annot$Tumor_Subtypes == "Lung"])

  p <- plot_heatmap(res_annot, "res_annot", subset_tumor = "Lung")
  expect_setequal(colnames(p@matrix), lung)
  expect_identical(names(p@top_annotation@anno_list), "Tumor_Subtypes")
  observed <- p@top_annotation@anno_list$Tumor_Subtypes@fun@var_env$value
  expect_equal(unname(observed), rep("Lung", length(lung)))
})

test_that(desc = "subset_tumor takes several values and shows 'annotation'", {
  data(res_annot)
  samp_annot(res_annot, "test") <- rep(c("A", "B", "C"), 4)
  annot <- samp_annot(res_annot)
  keep <- annot$Tumor_Subtypes == "Skin" | annot$test == "A"
  samples <- as.character(annot$Samples[keep])

  p <- plot_heatmap(res_annot, "res_annot", subset_tumor = c("Skin", "A"),
                    annotation = "test")
  expect_setequal(colnames(p@matrix), samples)
  expect_identical(names(p@top_annotation@anno_list), "test")
  observed <- p@top_annotation@anno_list$test@fun@var_env$value
  expected <- annot$test[match(colnames(p@matrix), annot$Samples)]
  expect_equal(unname(observed), expected)
})

test_that(desc = "subset_tumor names values not found in the annotations", {
  data(res_annot)
  expect_error(plot_heatmap(res_annot, "res_annot",
                            subset_tumor = c("Lung", "Brain")),
               regexp = "Brain")
})

test_that(desc = "scale = TRUE z-scores each signature across samples", {
  data(res_annot)
  p <- plot_heatmap(res_annot, "res_annot", scale = TRUE)
  observed <- p@matrix
  expect_equal(unname(rowMeans(observed)), rep(0, nrow(observed)))
  expect_equal(unname(apply(observed, 1, stats::sd)), rep(1, nrow(observed)))
})

test_that(desc = "scale = TRUE gives 0 for a constant signature", {
  data(res_annot)
  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  exposure_matrix[1, ] <- 5
  exposures(res_annot, "result", "SBS96", "res_annot") <- exposure_matrix
  p <- plot_heatmap(res_annot, "res_annot", scale = TRUE)
  expect_equal(unname(p@matrix[1, ]), rep(0, ncol(exposure_matrix)))
})

test_that(desc = "Proportional heatmap names samples with zero exposure", {
  data(res_annot)
  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  exposure_matrix[, 1] <- 0
  exposures(res_annot, "result", "SBS96", "res_annot") <- exposure_matrix
  expect_error(plot_heatmap(res_annot, "res_annot", proportional = TRUE),
               regexp = colnames(exposure_matrix)[1])
})

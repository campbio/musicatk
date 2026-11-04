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
  samp_annot(res_annot, "test") <- c("A", "B", "C", "A", "B", "C", "C")
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

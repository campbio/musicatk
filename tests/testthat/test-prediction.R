library(musicatk)

test_that("combine_predict_grid() keeps sample order and model fields", {
  data(musica_annot)
  data(cosmic_v2_sigs)
  set.seed(1)
  grid <- auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda",
    "Tumor_Subtypes",
    combine_res = FALSE, make_copy = TRUE
  )
  combined <- combine_predict_grid(grid, "SBS96", cosmic_v2_sigs,
    model_rename = "combined", make_copy = TRUE
  )
  model <- get_model(combined, "result", "SBS96", "combined")
  expect_identical(colnames(exposures(model)), sample_names(musica_annot))
  expect_identical(rownames(exposures(model)), colnames(signatures(model)))
  expect_equal(num_signatures(model), ncol(signatures(model)))
  expect_identical(model_id(model), "combined")
})

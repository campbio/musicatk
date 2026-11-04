library(musicatk)

test_that("combine_predict_grid() keeps sample order and model fields", {
  data(res_annot)
  model <- get_model(res_annot, "result", "SBS96", "res_annot")
  samples <- sample_names(res_annot)
  # Per-group models for interleaved groups, so binding by group is not
  # sample order
  group_a <- samples[c(1, 3, 5, 7)]
  group_b <- samples[c(2, 4, 6)]
  grid <- add_result(signatures(model), exposures(model)[, group_a],
    res_annot, "grid", "SBS96", "a",
    make_copy = TRUE
  )
  grid <- add_result(signatures(model), exposures(model)[, group_b],
    grid, "grid", "SBS96", "b",
    make_copy = TRUE
  )

  combined <- combine_predict_grid(grid, "SBS96", model,
    result_name = "grid", model_rename = "combined", make_copy = TRUE
  )
  combined_model <- get_model(combined, "grid", "SBS96", "combined")
  expect_identical(colnames(exposures(combined_model)), samples)
  expect_equal(
    exposures(combined_model),
    exposures(model)[sort(rownames(exposures(model))), samples]
  )
  expect_equal(
    num_signatures(combined_model),
    ncol(signatures(combined_model))
  )
  expect_identical(model_id(combined_model), "combined")
})

test_that("combine_predict_grid() stops when models share samples", {
  data(res_annot)
  model <- get_model(res_annot, "result", "SBS96", "res_annot")
  samples <- sample_names(res_annot)
  grid <- add_result(signatures(model), exposures(model)[, samples[1:4]],
    res_annot, "grid", "SBS96", "a",
    make_copy = TRUE
  )
  grid <- add_result(signatures(model), exposures(model)[, samples[3:7]],
    grid, "grid", "SBS96", "b",
    make_copy = TRUE
  )
  expect_error(
    combine_predict_grid(grid, "SBS96", model,
      result_name = "grid", make_copy = TRUE
    ),
    regexp = "more than one model"
  )
})

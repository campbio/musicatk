# A model_id already used in a result and modality stops with an error,
# unless overwrite = TRUE, which replaces that model in place.

clash_msg <- "already exists in result 'result', modality 'SBS96'"

discover_tiny <- function(musica, model_id, ...) {
  discover_signatures(musica, "SBS96", num_signatures = 2, nstart = 1,
                      model_id = model_id, make_copy = TRUE, ...)
}

test_that("discover_signatures() stops on a model_id already used", {
  musica <- discover_tiny(musica_sbs96_tiny, "m1")
  expect_error(discover_tiny(musica, "m1"), clash_msg)
  # The default ID clashes too
  musica <- discover_tiny(musica, NULL)
  expect_error(discover_tiny(musica, NULL), "'lda2' already exists")
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("m1", "lda2"))
})

test_that("discover_signatures() checks the model_id before fitting", {
  musica <- discover_tiny(musica_sbs96_tiny, "m1")
  local_mocked_bindings(LDA = function(...) stop("the fit ran"),
                        .package = "topicmodels")
  expect_error(discover_tiny(musica, "m1"), clash_msg)
  # The mock is in place: a new ID reaches the fit
  expect_error(discover_tiny(musica, "m2"), "the fit ran")
})

test_that("discover_signatures() with overwrite = TRUE replaces the model", {
  musica <- discover_tiny(musica_sbs96_tiny, "m1")
  musica <- discover_tiny(musica, "m2")
  musica <- discover_signatures(musica, "SBS96", num_signatures = 3,
                                nstart = 1, model_id = "m1",
                                make_copy = TRUE, overwrite = TRUE)
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("m1", "m2"))
  expect_identical(nrow(exposures(musica, "result", "SBS96", "m1")), 3L)
})

test_that("discover_signatures() still adds a new model_id", {
  musica <- discover_tiny(musica_sbs96_tiny, "m1")
  musica <- discover_tiny(musica, "m2")
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("m1", "m2"))
})

test_that("discover_signatures() in place stops on a clash", {
  musica <- musica_sbs96_tiny
  discover_signatures(musica, "SBS96", num_signatures = 2, nstart = 1,
                      model_id = "m1")
  expect_error(
    discover_signatures(musica, "SBS96", num_signatures = 2, nstart = 1,
                        model_id = "m1"),
    clash_msg
  )
  expect_identical(names(get_modality(musica, "result", "SBS96")), "m1")
})

predict_tiny <- function(musica, model_id, sigs = 1:3, ...) {
  predict_exposure(musica, "SBS96", cosmic_v2_sigs, algorithm = "lda",
                   model_id = model_id, signatures_to_use = sigs,
                   make_copy = TRUE, ...)
}

test_that("predict_exposure() stops on a model_id already used", {
  musica <- predict_tiny(musica_sbs96_tiny, "p1")
  expect_error(predict_tiny(musica, "p1"), clash_msg)
  musica <- predict_tiny(musica, NULL)
  expect_error(predict_tiny(musica, NULL), "'lda3_exp_pred' already exists")
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("p1", "lda3_exp_pred"))
})

test_that("predict_exposure() with overwrite = TRUE replaces the model", {
  musica <- predict_tiny(musica_sbs96_tiny, "p1")
  musica <- predict_tiny(musica, "p2")
  musica <- predict_tiny(musica, "p1", sigs = 1:5, overwrite = TRUE)
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("p1", "p2"))
  expect_identical(nrow(exposures(musica, "result", "SBS96", "p1")), 5L)
})

test_that("predict_exposure() still adds a new model_id", {
  musica <- predict_tiny(musica_sbs96_tiny, "p1")
  musica <- predict_tiny(musica, "p2")
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("p1", "p2"))
})

add_tiny <- function(musica, model_id, n_sigs = 2, ...) {
  sigs <- signatures(cosmic_v2_sigs)[, seq_len(n_sigs), drop = FALSE]
  exps <- matrix(1, nrow = n_sigs, ncol = 2,
                 dimnames = list(colnames(sigs), sample_names(musica)))
  add_result(sigs, exps, musica, result_name = "result", modality = "SBS96",
             model_id = model_id, make_copy = TRUE, ...)
}

test_that("add_result() stops on a model_id already used", {
  musica <- add_tiny(musica_sbs96_tiny, "a1")
  expect_error(add_tiny(musica, "a1", n_sigs = 4), clash_msg)
  expect_identical(nrow(exposures(musica, "result", "SBS96", "a1")), 2L)
})

test_that("add_result() with overwrite = TRUE replaces the model", {
  musica <- add_tiny(musica_sbs96_tiny, "a1")
  musica <- add_tiny(musica, "a2")
  musica <- add_tiny(musica, "a1", n_sigs = 4, overwrite = TRUE)
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("a1", "a2"))
  expect_identical(nrow(exposures(musica, "result", "SBS96", "a1")), 4L)
})

test_that("add_result() still adds a new model_id", {
  musica <- add_tiny(musica_sbs96_tiny, "a1")
  musica <- add_tiny(musica, "a2")
  expect_identical(names(get_modality(musica, "result", "SBS96")),
                   c("a1", "a2"))
})

test_that("overwrite must be a single TRUE or FALSE", {
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1)) {
    expect_error(add_tiny(musica_sbs96_tiny, "a1", overwrite = bad),
                 "'overwrite' must be TRUE or FALSE")
    expect_error(discover_tiny(musica_sbs96_tiny, "m1", overwrite = bad),
                 "'overwrite' must be TRUE or FALSE")
    expect_error(predict_tiny(musica_sbs96_tiny, "p1", overwrite = bad),
                 "'overwrite' must be TRUE or FALSE")
  }
})

test_that("auto_predict_grid() still adds .1 to a model_id already used", {
  musica <- auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda",
                              make_copy = TRUE)
  musica <- auto_predict_grid(musica, "SBS96", cosmic_v2_sigs, "lda",
                              make_copy = TRUE)
  ids <- names(get_modality(musica, "result", "SBS96"))
  expect_length(ids, 2)
  expect_identical(ids[2], paste0(ids[1], ".1"))
})

test_that("model_id must be a single string", {
  for (bad in list(NA_character_, c("a", "b"), character(0))) {
    expect_error(add_tiny(musica_sbs96_tiny, bad),
                 "'model_id' must be a single string")
    expect_error(predict_tiny(musica_sbs96_tiny, bad),
                 "'model_id' must be a single string")
  }
})

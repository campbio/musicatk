library(musicatk)

test_that("hyperparameter<- assigns on musica and result_collection", {
  data(res)
  new_hyper <- list(alpha = 0.5)
  hyperparameter(res, "result") <- new_hyper
  expect_identical(hyperparameter(res, "result"), new_hyper)

  collection <- result_list(res)[["result"]]
  hyperparameter(collection) <- list(beta = 2)
  expect_identical(hyperparameter(collection), list(beta = 2))
})

test_that("result_model setters work on musica, collection and model", {
  data(res)
  params <- S4Vectors::SimpleList(a = 1)
  intervals <- S4Vectors::SimpleList(b = 2)

  model <- get_model(res, "result", "SBS96", "res")
  num_signatures(model) <- 5
  other_parameters(model) <- params
  credible_intervals(model) <- intervals
  modality(model) <- "DBS78"
  model_id(model) <- "new_id"
  expect_equal(num_signatures(model), 5)
  expect_identical(other_parameters(model), params)
  expect_identical(credible_intervals(model), intervals)
  expect_identical(modality(model), "DBS78")
  expect_identical(model_id(model), "new_id")
  expect_error(model_id(model) <- c("a", "b"), regexp = "single value")

  collection <- result_list(res)[["result"]]
  num_signatures(collection, "SBS96", "res") <- 5
  other_parameters(collection, "SBS96", "res") <- params
  credible_intervals(collection, "SBS96", "res") <- intervals
  modality(collection, "SBS96", "res") <- "SBS96"
  expect_equal(num_signatures(collection, "SBS96", "res"), 5)
  expect_identical(other_parameters(collection, "SBS96", "res"), params)
  expect_identical(credible_intervals(collection, "SBS96", "res"), intervals)
  expect_identical(modality(collection, "SBS96", "res"), "SBS96")

  num_signatures(res, "result", "SBS96", "res") <- 5
  other_parameters(res, "result", "SBS96", "res") <- params
  credible_intervals(res, "result", "SBS96", "res") <- intervals
  modality(res, "result", "SBS96", "res") <- "SBS96"
  expect_equal(num_signatures(res, "result", "SBS96", "res"), 5)
  expect_identical(other_parameters(res, "result", "SBS96", "res"), params)
  expect_identical(credible_intervals(res, "result", "SBS96", "res"),
                   intervals)
  expect_identical(modality(res, "result", "SBS96", "res"), "SBS96")

  # The stored model is validated, not only the object that holds it
  expect_error(
    modality(collection, "SBS96", "res") <- c("A", "B"),
    regexp = "single value"
  )
  expect_error(
    modality(res, "result", "SBS96", "res") <- c("A", "B"),
    regexp = "single value"
  )
  expect_error(
    num_signatures(res, "result", "SBS96", "missing") <- 5,
    regexp = "not found"
  )
})

test_that("model_id<- renames the stored model so get_model() finds it", {
  data(res)
  collection <- result_list(res)[["result"]]
  model_id(collection, "SBS96", "res") <- "renamed"
  expect_identical(names(get_modality(collection, "SBS96")), "renamed")
  expect_identical(model_id(collection, "SBS96", "renamed"), "renamed")

  model_id(res, "result", "SBS96", "res") <- "renamed"
  model <- get_model(res, "result", "SBS96", "renamed")
  expect_identical(model_id(model), "renamed")
  expect_false("res" %in% names(get_modality(res, "result", "SBS96")))

  # Keeping the same id is allowed; taking another model's id is not
  model_id(res, "result", "SBS96", "renamed") <- "renamed"
  expect_identical(model_id(res, "result", "SBS96", "renamed"), "renamed")
  res <- add_result(signatures(model), exposures(model), res, "result",
                    "SBS96", "other", make_copy = TRUE)
  expect_error(
    model_id(res, "result", "SBS96", "other") <- "renamed",
    regexp = "already exists"
  )
  expect_error(
    model_id(res, "result", "SBS96", "missing") <- "x",
    regexp = "not found"
  )
})

test_that("modality<- moves the stored model to the new modality", {
  data(res)
  collection <- result_list(res)[["result"]]
  modality(collection, "SBS96", "res") <- "DBS78"
  expect_identical(modality(collection, "DBS78", "res"), "DBS78")
  expect_null(get_modality(collection, "SBS96"))

  modality(res, "result", "SBS96", "res") <- "DBS78"
  model <- get_model(res, "result", "DBS78", "res")
  expect_identical(modality(model), "DBS78")
  expect_null(get_modality(res, "result", "SBS96"))

  # Moving onto a model with the same id in the new modality is an error
  res <- add_result(signatures(model), exposures(model), res, "result",
                    "SBS96", "res", make_copy = TRUE)
  expect_error(
    modality(res, "result", "SBS96", "res") <- "DBS78",
    regexp = "already exists"
  )
})

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

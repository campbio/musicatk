library(musicatk)

# Example SBS96 table with the first sample's counts set to zero
.zero_sample_musica <- function() {
  counts <- get_count_table(
    extract_count_tables(musicatk::musica_sbs96)$SBS96
  )
  counts[, 1] <- 0
  create_musica_from_counts(counts, "SBS96")
}

test_that(".cosine() keeps a matrix for one signature (#162)", {
  sigs <- signatures(musicatk::cosmic_v2_sigs)
  res <- musicatk:::.cosine(sigs[, 1, drop = FALSE], sigs)
  expect_identical(dim(res), c(1L, ncol(sigs)))

  one_sig <- predict_exposure(musicatk::musica_sbs96, "SBS96",
    musicatk::cosmic_v2_sigs, "decompTumor2Sig",
    signatures_to_use = 1, model_id = "one", make_copy = TRUE
  )
  comparison <- compare_cosmic_v2(one_sig, "one", threshold = 0.99)
  expect_identical(comparison$y_sig_name, "SBS1")
})

test_that("model_id de-duplication counts past 10 (#148)", {
  m <- musicatk::musica_sbs96
  for (i in seq_len(12)) {
    m <- suppressMessages(predict_exposure(m, "SBS96",
      musicatk::cosmic_v2_sigs, "decompTumor2Sig",
      signatures_to_use = 1:2, model_id = "p", make_copy = TRUE
    ))
  }
  expect_identical(
    names(get_modality(m, "result", "SBS96")),
    c("p", paste0("p.", 1:11))
  )
})

test_that("predict_exposure() matches motifs by name (#142)", {
  sigs <- musicatk::cosmic_v2_sigs
  permuted <- sigs
  signatures(permuted) <- signatures(sigs)[rev(seq_len(96)), ]
  expected <- predict_exposure(musicatk::musica_sbs96, "SBS96", sigs,
    "decompTumor2Sig",
    signatures_to_use = 1:5, model_id = "p", make_copy = TRUE
  )
  got <- predict_exposure(musicatk::musica_sbs96, "SBS96", permuted,
    "decompTumor2Sig",
    signatures_to_use = 1:5, model_id = "p", make_copy = TRUE
  )
  expect_equal(
    exposures(got, "result", "SBS96", "p"),
    exposures(expected, "result", "SBS96", "p")
  )
  expect_equal(
    signatures(got, "result", "SBS96", "p"),
    signatures(expected, "result", "SBS96", "p")
  )

  expect_error(
    predict_exposure(musicatk::musica_sbs96, "SBS96",
      musicatk::cosmic_v3_dbs_sigs, "lda",
      make_copy = TRUE
    ),
    regexp = "motifs"
  )
  expect_error(
    predict_exposure(musicatk::dbs_musica, "DBS78",
      musicatk::cosmic_v3_dbs_sigs, "decompTumor2Sig",
      make_copy = TRUE
    ),
    regexp = "SBS96"
  )
})

test_that("zero-count samples get zero exposures (#143)", {
  zero <- .zero_sample_musica()
  samples <- sample_names(zero)

  for (algorithm in c("lda", "nmf")) {
    expect_message(
      m <- discover_signatures(zero, "SBS96", 2, algorithm,
        seed = 1, nstart = 1, make_copy = TRUE
      ),
      regexp = "no mutations"
    )
    e <- exposures(m, "result", "SBS96", paste0(algorithm, 2))
    expect_identical(colnames(e), samples)
    expect_true(all(e[, 1] == 0))
    expect_true(all(colSums(e)[-1] > 0))
  }

  for (algorithm in c("lda", "decompTumor2Sig")) {
    expect_message(
      m <- predict_exposure(zero, "SBS96", musicatk::cosmic_v2_sigs,
        algorithm,
        signatures_to_use = 1:3, model_id = "p", make_copy = TRUE
      ),
      regexp = "no mutations"
    )
    e <- exposures(m, "result", "SBS96", "p")
    expect_identical(colnames(e), samples)
    expect_true(all(e[, 1] == 0))
  }

  grid <- suppressMessages(generate_result_grid(zero, "SBS96", "lda",
    k_start = 2, k_end = 2, seed = 1, make_copy = TRUE
  ))
  model <- get_model(grid, "result_grid", "SBS96", "lda2")
  expect_identical(colnames(exposures(model)), samples)
  expect_true(all(exposures(model)[, 1] == 0))
  expect_true(is.finite(metrics(model)$reconstruction_error))

  # compare_k_vals() prints a plot; don't write Rplots.pdf
  withr::local_pdf(NULL)
  withr::local_seed(1)
  k_vals <- suppressMessages(compare_k_vals(zero, "SBS96",
    reps = 1, min_k = 2, max_k = 3
  ))
  expect_true(all(is.finite(k_vals$error)))

  counts <- get_count_table(extract_count_tables(zero)$SBS96)
  counts[, -2] <- 0
  one_left <- create_musica_from_counts(counts, "SBS96")
  expect_error(
    discover_signatures(one_left, "SBS96", 2, seed = 1, nstart = 1,
      make_copy = TRUE
    ),
    regexp = "at least two samples with mutations"
  )
})

test_that("generate_result_grid() splits by the annotation argument (#140)", {
  grid <- generate_result_grid(musicatk::musica_annot, "SBS96", "lda",
    annotation = "Tumor_Subtypes", k_start = 2, k_end = 2, seed = 1,
    make_copy = TRUE
  )
  models <- get_modality(grid, "result_grid", "SBS96")
  expect_identical(sort(names(models)), c("lda2_Lung", "lda2_Skin"))
  annot <- samp_annot(musicatk::musica_annot)
  lung <- annot$Samples[annot$Tumor_Subtypes == "Lung"]
  skin <- annot$Samples[annot$Tumor_Subtypes == "Skin"]
  expect_setequal(colnames(exposures(models[["lda2_Lung"]])), lung)
  expect_setequal(colnames(exposures(models[["lda2_Skin"]])), skin)

  expect_error(
    generate_result_grid(musicatk::musica_annot, "SBS96", "lda",
      annotation = "Not_A_Column", k_start = 2, k_end = 2,
      make_copy = TRUE
    ),
    regexp = "Not_A_Column"
  )
})

test_that("generate_result_grid() re-run stores the new fit (#141)", {
  grid <- generate_result_grid(musicatk::musica_sbs96, "SBS96", "lda",
    k_start = 2, k_end = 2, seed = 1, make_copy = TRUE
  )
  first <- exposures(grid, "result_grid", "SBS96", "lda2")
  expect_message(
    grid <- generate_result_grid(grid, "SBS96", "lda",
      k_start = 2, k_end = 2, seed = 2, make_copy = TRUE
    ),
    regexp = "lda2.1"
  )
  models <- get_modality(grid, "result_grid", "SBS96")
  expect_setequal(names(models), c("lda2", "lda2.1"))
  expect_identical(exposures(models[["lda2"]]), first)
  expect_false(isTRUE(all.equal(exposures(models[["lda2.1"]]), first)))
  expect_true(is.finite(metrics(models[["lda2.1"]])$reconstruction_error))
})

test_that("auto_predict_grid() works with any result_name (#145)", {
  m <- auto_predict_grid(musicatk::musica_annot, "SBS96",
    musicatk::cosmic_v2_sigs, "decompTumor2Sig",
    result_name = "pred", make_copy = TRUE
  )
  expect_length(get_modality(m, "pred", "SBS96"), 1)
})

test_that("auto_predict_grid() accepts a factor annotation (#146)", {
  from_factor <- musicatk::musica_annot
  samp_annot(from_factor, "Tumor_Subtypes") <-
    factor(samp_annot(from_factor)$Tumor_Subtypes)
  from_factor <- auto_predict_grid(from_factor, "SBS96",
    musicatk::cosmic_v2_sigs, "decompTumor2Sig",
    sample_annotation = "Tumor_Subtypes", verbose = FALSE, make_copy = TRUE
  )
  from_character <- auto_predict_grid(musicatk::musica_annot, "SBS96",
    musicatk::cosmic_v2_sigs, "decompTumor2Sig",
    sample_annotation = "Tumor_Subtypes", verbose = FALSE, make_copy = TRUE
  )
  expect_identical(
    get_modality(from_factor, "result", "SBS96"),
    get_modality(from_character, "result", "SBS96")
  )
})

test_that("auto_subset_sigs() explains when no signature passes (#147)", {
  expect_error(
    musicatk:::auto_subset_sigs(musicatk::musica_annot, "SBS96",
      musicatk::cosmic_v2_sigs, "decompTumor2Sig",
      min_exists = 2, proportion_samples = 2, rare_exposure = 2
    ),
    regexp = "No signature passed"
  )
})

test_that("compare_k_vals() checks error_type before fitting (#149)", {
  messages <- character(0)
  expect_error(
    withCallingHandlers(
      compare_k_vals(musicatk::musica_sbs96, "SBS96",
        reps = 1, min_k = 2, max_k = 2, error_type = "bogus"
      ),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    ),
    regexp = "error_type"
  )
  # no fit was started
  expect_length(messages, 0)
})

test_that("LDA prediction is seeded and keeps the caller's RNG (#152)", {
  withr::local_seed(42)
  before <- get(".Random.seed", envir = globalenv())
  first <- predict_exposure(musicatk::musica_sbs96, "SBS96",
    musicatk::cosmic_v2_sigs, "lda",
    signatures_to_use = 1:5, model_id = "p", make_copy = TRUE
  )
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
  second <- predict_exposure(musicatk::musica_sbs96, "SBS96",
    musicatk::cosmic_v2_sigs, "lda",
    signatures_to_use = 1:5, model_id = "p", make_copy = TRUE
  )
  expect_identical(
    exposures(first, "result", "SBS96", "p"),
    exposures(second, "result", "SBS96", "p")
  )
})

test_that("LDA prediction of one sample doesn't warn (#152)", {
  counts <- get_count_table(
    extract_count_tables(musicatk::musica_sbs96)$SBS96
  )
  one <- create_musica_from_counts(counts[, 1, drop = FALSE], "SBS96")
  expect_no_warning(predict_exposure(one, "SBS96",
    musicatk::cosmic_v2_sigs, "lda",
    signatures_to_use = 1:5, model_id = "p", make_copy = TRUE
  ))
})

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

test_that("compare_k_vals() checks error_type before fitting (#149)", {
  expect_error(
    compare_k_vals(musicatk::musica_sbs96, "SBS96",
      reps = 1, min_k = 2, max_k = 2, error_type = "bogus"
    ),
    regexp = "error_type"
  )
})

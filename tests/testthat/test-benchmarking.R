library(musicatk)

true_sig_names <- c(
  "SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26"
)

test_that("create_benchmark() stores true exposures matched by sample name", {
  data(synthetic_breast_counts)
  data(synthetic_breast_true_exposures)
  data(cosmic_v2_sigs)
  true_signatures <- signatures(cosmic_v2_sigs)[, true_sig_names]
  expected <- t(as.matrix(synthetic_breast_true_exposures))
  colnames(expected) <- rownames(synthetic_breast_true_exposures)

  # Rows (samples) in a different order from the count table
  shuffled <- synthetic_breast_true_exposures[
    rev(seq_len(nrow(synthetic_breast_true_exposures))),
  ]
  bench <- create_benchmark(true_signatures, shuffled,
    synthetic_breast_counts, "SBS96"
  )
  truth_musica <- ground_truth(bench)
  truth <- exposures(truth_musica, "ground_truth", "SBS96", "ground_truth")
  expect_identical(colnames(truth), colnames(synthetic_breast_counts))
  expect_equal(truth, expected)

  # Signatures x samples input is accepted too
  bench_t <- create_benchmark(true_signatures, expected,
    synthetic_breast_counts, "SBS96"
  )
  truth_musica_t <- ground_truth(bench_t)
  truth_t <- exposures(truth_musica_t, "ground_truth", "SBS96", "ground_truth")
  expect_equal(truth_t, expected)

  renamed <- synthetic_breast_true_exposures
  rownames(renamed)[1] <- "not_a_sample"
  expect_error(
    create_benchmark(
      true_signatures, renamed, synthetic_breast_counts, "SBS96"
    ),
    regexp = "sample names of 'true_exposures'"
  )
})

test_that("benchmark_plot_exposures() pairs exposures by sample name", {
  data(synthetic_breast_counts)
  data(synthetic_breast_true_exposures)
  data(example_predicted_sigs)
  data(example_predicted_exp)
  data(cosmic_v2_sigs)
  true_signatures <- signatures(cosmic_v2_sigs)[, true_sig_names]
  bench <- create_benchmark(true_signatures, synthetic_breast_true_exposures,
    synthetic_breast_counts, "SBS96"
  )
  prediction <- create_musica_from_counts(synthetic_breast_counts, "SBS96")
  add_result(
    example_predicted_sigs, example_predicted_exp, prediction,
    "result", "SBS96", "example_k8"
  )
  bench <- suppressWarnings(benchmark(bench, prediction, "example_k8",
    "SBS96",
    method_id = "example", plot = FALSE, make_copy = TRUE
  ))
  p <- benchmark_plot_exposures(bench, "example", "Initial")
  comparison <- benchmark_get_comparison(bench, "example", "Initial")
  first <- p$data[p$data$Sig == comparison$y_sig_name[1], ]
  predicted_sig <- comparison[1, 4]
  pred_exp <- example_predicted_exp[predicted_sig, ]
  true_sig <- comparison$y_sig_name[1]
  true_exp <- synthetic_breast_true_exposures[names(pred_exp), true_sig]
  expect_equal(first$Predicted, unname(pred_exp))
  expect_equal(first$True, true_exp)
})

test_that("benchmark() handles one unchanged signature plus a duplicate", {
  data(cosmic_v2_sigs)
  set.seed(1)
  true_signatures <- signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")]
  samples <- paste0("sample", seq_len(20))
  true_exposures <- matrix(
    sample(100:1000, 40, replace = TRUE),
    nrow = 2, dimnames = list(colnames(true_signatures), samples)
  )
  counts <- round(true_signatures %*% true_exposures)
  bench <- create_benchmark(true_signatures, true_exposures, counts, "SBS96")

  # SBS1 predicted once; SBS2 split into two near-copies (a duplicate)
  noise <- matrix(runif(96 * 2, 0, 0.002), ncol = 2)
  split_sigs <- true_signatures[, c("SBS2", "SBS2")] + noise
  split_sigs <- sweep(split_sigs, 2, colSums(split_sigs), "/")
  predicted_sigs <- cbind(true_signatures[, "SBS1"], split_sigs)
  colnames(predicted_sigs) <- c("Signature1", "Signature2", "Signature3")
  predicted_exposures <- rbind(
    true_exposures["SBS1", ],
    true_exposures["SBS2", ] / 2,
    true_exposures["SBS2", ] / 2
  )
  rownames(predicted_exposures) <- colnames(predicted_sigs)
  prediction <- create_musica_from_counts(counts, "SBS96")
  add_result(predicted_sigs, predicted_exposures, prediction, "result",
    "SBS96", "split"
  )

  bench <- suppressWarnings(benchmark(bench, prediction, "split", "SBS96",
    method_id = "split", plot = FALSE, make_copy = TRUE
  ))
  intermediate <- intermediate_pred(benchmark_get_entry(bench, "split"))
  expect_identical(
    rownames(exposures(intermediate)),
    colnames(signatures(intermediate))
  )
  expect_true("Signature1" %in% colnames(signatures(intermediate)))
})

test_that("Benchmarks saved with older versions give a clear error", {
  data(full_benchmark_example)
  truth <- get_model(ground_truth(full_benchmark_example), "ground_truth",
    "SBS96", "ground_truth"
  )
  samples <- colnames(exposures(truth))
  # Older versions stored the truth samples x signatures, without names
  legacy <- truth
  legacy@exposures <- unname(t(exposures(truth)))
  colnames(legacy@exposures) <- rownames(exposures(truth))
  expect_error(
    .true_exposure(legacy, "SBS1", samples),
    regexp = "older version of musicatk"
  )
})

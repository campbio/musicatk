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

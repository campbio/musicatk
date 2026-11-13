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

test_that("Benchmark functions accept subclasses of full_benchmark (#94)", {
  data(full_benchmark_example)
  setClass("full_benchmark_sub", contains = "full_benchmark",
    where = globalenv()
  )
  withr::defer(removeClass("full_benchmark_sub", where = globalenv()))
  sub <- new("full_benchmark_sub", full_benchmark_example)
  entry <- benchmark_get_entry(sub, "example_k8")
  expect_s4_class(entry, "single_benchmark")
  expect_s3_class(
    benchmark_get_comparison(sub, "example_k8", "Final"), "data.frame"
  )
  expect_error(
    benchmark_get_entry(list(), "example_k8"),
    regexp = "is not a 'full_benchmark' object"
  )
})

test_that("Stage aliases work in every benchmark accessor (#193)", {
  data(full_benchmark_example)
  final <- benchmark_get_comparison(full_benchmark_example, "example_k8",
    "Final"
  )
  for (alias in c("final", "Fin", "fin")) {
    expect_identical(
      benchmark_get_comparison(full_benchmark_example, "example_k8", alias),
      final
    )
  }
  entry <- benchmark_get_entry(full_benchmark_example, "example_k8")
  expect_identical(
    benchmark_get_comparison(full_benchmark_example, "example_k8", "init"),
    initial_comparison(entry)
  )
  expect_identical(
    benchmark_get_comparison(full_benchmark_example, "example_k8", "inter"),
    intermediate_comparison(entry)
  )
  expect_error(
    benchmark_get_comparison(full_benchmark_example, "example_k8", "last"),
    regexp = "'prediction' must be one of"
  )
})

test_that("create_benchmark() stores the truth under the table name (#195)", {
  data(cosmic_v2_sigs)
  data(cosmic_v3_dbs_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")])
  setup <- setup_benchmark(truth, truth$signatures, truth$exposures,
    modality = "SBS96", variant_class = "SBS"
  )
  truth_musica <- ground_truth(setup$bench)
  expect_false(is.null(get_modality(truth_musica, "ground_truth", "SBS96")))
  bench <- run_benchmark(setup$bench, setup$prediction, "SBS96")
  expect_named(indv_benchmarks(bench), "m")
  expect_error(
    run_benchmark(setup$bench, setup$prediction, "DBS78"),
    regexp = "'prediction_modality' .*DBS78.* ground truth .*SBS96"
  )

  dbs_truth <- sim_truth(signatures(cosmic_v3_dbs_sigs)[, c("DBS2", "DBS4")])
  dbs_setup <- setup_benchmark(dbs_truth, dbs_truth$signatures,
    dbs_truth$exposures,
    modality = "DBS78"
  )
  dbs_bench <- run_benchmark(dbs_setup$bench, dbs_setup$prediction, "DBS78")
  expect_named(indv_benchmarks(dbs_bench), "m")

  expect_error(
    create_benchmark(truth$signatures, truth$exposures, truth$counts,
      "IND83"
    ),
    regexp = "'variant_class' must be"
  )
})

test_that("benchmark(plot = FALSE) draws nothing (#196)", {
  data(cosmic_v2_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")])
  setup <- setup_benchmark(truth, truth$signatures, truth$exposures)
  local_mocked_bindings(
    .plot_compare_result_signatures = function(...) stop("drew a plot")
  )
  expect_no_error(run_benchmark(setup$bench, setup$prediction))
})

test_that("Benchmark summaries are data.frames with typed columns (#201)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS2", "SBS3", "SBS8")])
  # P1 is SBS1 and P2 a composite of SBS3 and SBS8; SBS2 is missed
  pred <- withr::with_seed(5, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = noisy_sig(rowMeans(sigs[, c("SBS3", "SBS8")]))
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = colSums(truth$exposures[c("SBS3", "SBS8"), ])
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)

  # a global variable with this name must not change the counts
  assign("duplicated_signature_components", "P2", envir = globalenv())
  withr::defer(rm("duplicated_signature_components", envir = globalenv()))
  expect_no_warning(bench <- run_benchmark(setup$bench, setup$prediction))

  method_view <- method_view_summary(bench)
  expect_s3_class(method_view, "data.frame")
  expect_identical(rownames(method_view), "m")
  expect_type(method_view[["Num Found"]], "integer")
  expect_type(method_view[["Mean CS"]], "double")
  expect_type(method_view[["Sigs Found"]], "character")
  expect_identical(method_view[["Num Composites"]], 1L)
  expect_identical(method_view[["Num Duplicates"]], 0L)
  expect_identical(
    single_summary(benchmark_get_entry(bench, "m")), method_view
  )

  sig_view <- sig_view_summary(bench)
  expect_s3_class(sig_view, "data.frame")
  expect_setequal(rownames(sig_view), colnames(truth$signatures))
  expect_type(sig_view[["Times Found"]], "integer")
  expect_type(sig_view[["Max CS"]], "double")
  expect_type(sig_view[["Methods Found"]], "character")
  expect_identical(sig_view["SBS2", "Times Missed"], 1L)
  expect_true(is.na(sig_view["SBS2", "Min CS"]))

  # A full_benchmark from an older version, with matrix summaries
  legacy <- setup$bench
  attr(legacy, "method_view_summary") <- matrix(nrow = 0, ncol = 0)
  expect_error(
    run_benchmark(legacy, setup$prediction),
    regexp = "older version of musicatk"
  )
})

test_that("Initial RE follows re_method (#199)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS2")])
  pred <- withr::with_seed(6, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]), P2 = noisy_sig(sigs[, "SBS2"])
  ))
  pred_exposures <- truth$exposures
  rownames(pred_exposures) <- colnames(pred)
  setup <- setup_benchmark(truth, pred, pred_exposures)

  prop <- run_benchmark(setup$bench, setup$prediction)
  expected_prop <- sum(abs(
    prop.table(truth$counts, 2) - pred %*% prop.table(pred_exposures, 2)
  ))
  expect_equal(
    method_view_summary(prop)[["Initial RE"]], round(expected_prop, 3)
  )

  raw <- run_benchmark(setup$bench, setup$prediction, re_method = "raw")
  expected_raw <- sum(abs(truth$counts - pred %*% pred_exposures))
  expect_equal(
    method_view_summary(raw)[["Initial RE"]], round(expected_raw, 3)
  )
})

test_that("A prediction with no match is recorded as a total miss (#200)", {
  data(cosmic_v2_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2", "SBS3")])
  # one mutation type each, unlike any true signature
  pred <- matrix(0, 96, 3,
    dimnames = list(rownames(truth$signatures), c("P1", "P2", "P3"))
  )
  pred[1, "P1"] <- 1
  pred[50, "P2"] <- 1
  pred[90, "P3"] <- 1
  pred_exposures <- truth$exposures
  rownames(pred_exposures) <- colnames(pred)
  setup <- setup_benchmark(truth, pred, pred_exposures)

  expect_no_warning(bench <- run_benchmark(setup$bench, setup$prediction))
  summary <- method_view_summary(bench)
  expect_identical(summary[["Num Found"]], 0L)
  expect_identical(summary[["Num Missing"]], 3L)
  expect_identical(summary[["Num Spurious"]], 3L)
  expect_identical(summary[["Num Direct Matched"]], 0L)
  expect_true(is.na(summary[["Mean CS"]]))
  expect_true(is.na(summary[["Min CS"]]))
  expect_identical(summary[["Sigs Found"]], "")
  expect_equal(nrow(benchmark_get_comparison(bench, "m", "Final")), 0)
  expect_error(
    benchmark_plot_comparison(bench, "m", "Final"),
    regexp = "No predicted signature matches a true signature .*Final"
  )
  expect_error(
    benchmark_plot_exposures(bench, "m", "Initial"),
    regexp = "No predicted signature matches a true signature .*Initial"
  )
})

test_that("A predicted signature joins one duplicate group only (#188)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS3", "SBS8")])
  exposures <- truth$exposures
  # P1 and P2 duplicate SBS3; P4 matches SBS3 (best) and SBS8
  pred <- withr::with_seed(2, cbind(
    P0 = noisy_sig(sigs[, "SBS1"]),
    P1 = noisy_sig(sigs[, "SBS3"]),
    P2 = noisy_sig(sigs[, "SBS3"]),
    P3 = noisy_sig(sigs[, "SBS8"]),
    P4 = noisy_sig((sigs[, "SBS3"] + sigs[, "SBS8"]) / 2)
  ))
  pred_exposures <- rbind(
    P0 = exposures["SBS1", ],
    P1 = exposures["SBS3", ] * 0.4,
    P2 = exposures["SBS3", ] * 0.4,
    P3 = exposures["SBS8", ] * 0.8,
    P4 = (exposures["SBS3", ] + exposures["SBS8", ]) * 0.2
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)
  entry <- benchmark_get_entry(bench, "m")

  expect_equal(
    sum(exposures(intermediate_pred(entry))),
    sum(exposures(initial_pred(entry)))
  )
  intermediate_sigs <- colnames(signatures(intermediate_pred(entry)))
  expect_setequal(
    intermediate_sigs, c("P0", "P3", "Merged Signature (P1.P2.P4)")
  )
  summary <- method_view_summary(bench)
  expect_identical(summary[["Num Direct Matched"]], 2L)
  expect_identical(summary[["Num Duplicates"]], 2L)
  expect_identical(summary[["Num Dup/Comp"]], 1L)

  # ggplot2 warns that the fitted lines leave the plots' fixed axes
  plots <- suppressWarnings(benchmark_plot_duplicate_exposures(bench, "m"))
  expect_length(plots, 1)
})

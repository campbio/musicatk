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
  for (alias in c("final", "Fin", "fin", "FINAL")) {
    expect_identical(
      benchmark_get_comparison(full_benchmark_example, "example_k8", alias),
      final
    )
  }
  entry <- benchmark_get_entry(full_benchmark_example, "example_k8")
  expect_identical(
    benchmark_get_comparison(full_benchmark_example, "example_k8", "Init"),
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
  expect_identical(rownames(sig_view), colnames(truth$signatures))
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
  # SBS3 is matched by P1, P2 and P4, and SBS8 by P3 and P4: P1, P2 and P3
  # are duplicates, P4 is both a duplicate and a composite, P0 is direct
  summary <- method_view_summary(bench)
  expect_identical(summary[["Num Direct Matched"]], 1L)
  expect_identical(summary[["Num Duplicates"]], 3L)
  expect_identical(summary[["Num Dup/Comp"]], 1L)

  # ggplot2 warns that the fitted lines leave the plots' fixed axes
  plots <- suppressWarnings(benchmark_plot_duplicate_exposures(bench, "m"))
  expect_length(plots, 1)
})

# Benchmark of SBS1, SBS2 (predicted as such) and a composite 'P3' of the
# true signatures 'parts'
composite_benchmark <- function(parts, sd = 0.0005, seed = 3) {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS2", parts)])
  pred <- withr::with_seed(seed, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = noisy_sig(sigs[, "SBS2"]),
    P3 = noisy_sig(rowMeans(sigs[, parts]), sd)
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = truth$exposures["SBS2", ],
    P3 = colSums(truth$exposures[parts, ])
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  return(run_benchmark(setup$bench, setup$prediction))
}

test_that("A composite of three true signatures is split in three (#189)", {
  bench <- composite_benchmark(c("SBS3", "SBS5", "SBS8"))
  entry <- benchmark_get_entry(bench, "m")
  final <- final_pred(entry)
  pieces <- c("SBS3_like", "SBS5_like", "SBS8_like")
  expect_setequal(colnames(signatures(final)), c("P1", "P2", pieces))
  expect_equal(
    colSums(exposures(final)[pieces, ]),
    exposures(intermediate_pred(entry))["P3", ]
  )
})

test_that("Composites of DBS78 signatures are split (#190)", {
  data(cosmic_v3_dbs_sigs)
  dbs <- signatures(cosmic_v3_dbs_sigs)
  truth <- sim_truth(dbs[, c("DBS2", "DBS3", "DBS10")])
  pred <- withr::with_seed(7, cbind(
    P1 = noisy_sig(dbs[, "DBS2"]),
    P2 = noisy_sig(0.6 * dbs[, "DBS3"] + 0.4 * dbs[, "DBS10"], 0.0005)
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["DBS2", ],
    P2 = colSums(truth$exposures[c("DBS3", "DBS10"), ])
  )
  setup <- setup_benchmark(truth, pred, pred_exposures, modality = "DBS78")
  bench <- run_benchmark(setup$bench, setup$prediction, "DBS78")
  final <- final_pred(benchmark_get_entry(bench, "m"))
  expect_setequal(
    colnames(signatures(final)), c("P1", "DBS3_like", "DBS10_like")
  )
})

test_that("Composite pieces keep the composite's noise (#191)", {
  bench <- composite_benchmark(c("SBS3", "SBS8"), sd = 0.004, seed = 5)
  entry <- benchmark_get_entry(bench, "m")
  final_comp <- benchmark_get_comparison(bench, "m", "Final")
  pieces <- final_comp[grepl("_like$", final_comp$x_sig_name), ]
  expect_setequal(pieces$x_sig_name, c("SBS3_like", "SBS8_like"))
  # pieces are not copies of the true signatures
  expect_true(all(pieces$cosine < 0.999))

  # splitting doesn't change the reconstruction
  reconstruct <- function(model) {
    signatures(model) %*%
      exposures(model)[colnames(signatures(model)), , drop = FALSE]
  }
  expect_equal(
    reconstruct(final_pred(entry)), reconstruct(intermediate_pred(entry))
  )
  # each piece gets a fixed fraction of the composite's exposure
  final_exposures <- exposures(final_pred(entry))
  composite <- exposures(intermediate_pred(entry))["P3", ]
  fraction <- final_exposures["SBS3_like", ] / composite
  expect_equal(unname(fraction), rep(fraction[[1]], length(fraction)))
  summary <- method_view_summary(bench)
  expect_equal(summary[["Initial RE"]], summary[["Final RE"]])
  best <- tapply(final_comp$cosine, final_comp$y_sig_name, max)
  expect_equal(summary[["Mean CS"]], round(mean(best), 3))

  # ggplot2 warns that the fitted lines leave the plots' fixed axes
  plots <- suppressWarnings(benchmark_plot_composite_exposures(bench, "m"))
  expect_length(plots, 1)
})

test_that("Splitting an unnormalised composite keeps the reconstruction", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS3", "SBS8")])
  # P2's column sums to 2, as raw NMF signatures may
  pred <- withr::with_seed(8, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = 2 * noisy_sig(rowMeans(sigs[, c("SBS3", "SBS8")]), 0.0005)
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = colSums(truth$exposures[c("SBS3", "SBS8"), ]) / 2
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction, re_method = "raw")
  entry <- benchmark_get_entry(bench, "m")
  reconstruct <- function(model) {
    signatures(model) %*%
      exposures(model)[colnames(signatures(model)), , drop = FALSE]
  }
  expect_equal(
    reconstruct(final_pred(entry)), reconstruct(intermediate_pred(entry))
  )
})

test_that("A composite leaving one signature unchanged keeps its name (#192)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS3", "SBS8")])
  pred <- withr::with_seed(8, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = noisy_sig(rowMeans(sigs[, c("SBS3", "SBS8")]), 0.0005)
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = colSums(truth$exposures[c("SBS3", "SBS8"), ])
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)
  final <- final_pred(benchmark_get_entry(bench, "m"))
  expect_setequal(
    colnames(signatures(final)), c("P1", "SBS3_like", "SBS8_like")
  )
  expect_identical(rownames(exposures(final)), colnames(signatures(final)))
})

test_that("Summaries count each true signature once, in Final (#198)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS3", "SBS8")])
  # P3, a composite of SBS3 and SBS8, is split; on devel, its SBS3 piece and
  # P2 then both matched SBS3 in the Final stage and SBS3 counted twice
  pred <- withr::with_seed(4, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = noisy_sig(sigs[, "SBS3"]),
    P3 = noisy_sig(0.45 * sigs[, "SBS3"] + 0.55 * sigs[, "SBS8"], 0.001)
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = truth$exposures["SBS3", ] * 0.5,
    P3 = truth$exposures["SBS3", ] * 0.5 + truth$exposures["SBS8", ]
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)
  final <- benchmark_get_comparison(bench, "m", "Final")
  best <- tapply(final$cosine, final$y_sig_name, max)
  summary <- method_view_summary(bench)
  expect_identical(summary[["Num Found"]], 3L)
  expect_identical(summary[["Num Missing"]], 0L)
  expect_identical(summary[["Sigs Found"]], "SBS1, SBS3, SBS8")
  expect_equal(summary[["Mean CS"]], round(mean(best), 3))
  expect_equal(summary[["Min CS"]], round(min(best), 3))

  sig_view <- sig_view_summary(bench)
  expect_identical(sig_view["SBS3", "Times Found"], 1L)
  expect_equal(sig_view["SBS3", "Mean CS"], round(best[["SBS3"]], 3))

  # a true signature matched twice counts once, at its best cosine
  doubled <- rbind(final, final[final$y_sig_name == "SBS3", ])
  doubled$cosine[nrow(doubled)] <- 0.85
  expect_equal(
    .best_cosines(doubled, c("SBS1", "SBS3", "SBS8", "SBS13")),
    c(SBS1 = best[["SBS1"]], SBS3 = best[["SBS3"]], SBS8 = best[["SBS8"]])
  )
})

test_that("predict_and_benchmark() benchmarks the model it discovers (#197)", {
  data(cosmic_v2_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")])
  bench <- create_benchmark(truth$signatures, truth$exposures, truth$counts,
    "SBS96"
  )
  musica <- create_musica_from_counts(truth$counts, "SBS96")
  # an earlier lda2, made with another seed
  suppressMessages(discover_signatures(musica, "SBS96", 2, "lda",
    seed = 99, nstart = 1
  ))
  expected <- suppressMessages(discover_signatures(
    create_musica_from_counts(truth$counts, "SBS96"), "SBS96", 2, "lda",
    seed = 1, nstart = 1, make_copy = TRUE
  ))

  utils::capture.output(suppressMessages(
    result <- predict_and_benchmark(musica, "SBS96", 2, 2, bench, "lda",
      nstart = 1
    )
  ))
  expect_named(indv_benchmarks(result), "lda2.1")
  initial <- initial_pred(benchmark_get_entry(result, "lda2.1"))
  expect_identical(model_id(initial), "lda2.1")
  expect_equal(
    signatures(initial), signatures(expected, "result", "SBS96", "lda2")
  )
})

test_that("benchmark() makes method_id unique past eleven repeats (#148)", {
  data(cosmic_v2_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")])
  setup <- setup_benchmark(truth, truth$signatures, truth$exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)
  entry <- benchmark_get_entry(bench, "m")
  repeats <- rep(list(entry), 11)
  names(repeats) <- c("m", paste0("m.", 1:10))
  indv_benchmarks(bench) <- repeats
  bench <- run_benchmark(bench, setup$prediction)
  expect_identical(utils::tail(names(indv_benchmarks(bench)), 1), "m.11")
})

test_that("benchmark_plot_comparison() plots a stage with matches", {
  data(full_benchmark_example)
  path <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer(grDevices::dev.off())
  expect_no_error(
    benchmark_plot_comparison(full_benchmark_example, "example_k8", "final")
  )
})

test_that("predict_and_benchmark() checks the benchmark before discovery", {
  data(cosmic_v2_sigs)
  truth <- sim_truth(signatures(cosmic_v2_sigs)[, c("SBS1", "SBS2")])
  bench <- create_benchmark(truth$signatures, truth$exposures, truth$counts,
    "SBS96"
  )
  musica <- create_musica_from_counts(truth$counts, "SBS96")
  local_mocked_bindings(
    discover_signatures = function(...) stop("ran discovery")
  )
  expect_error(
    predict_and_benchmark(musica, "DBS78", 2, 2, bench),
    regexp = "'prediction_modality' is 'DBS78'"
  )
  expect_error(
    predict_and_benchmark(musica, "SBS96", 2, 2, list()),
    regexp = "'full_benchmark' must be a 'full_benchmark' object"
  )
})

test_that("A composite piece and a duplicate prediction are merged in Final", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS3", "SBS8")])
  # P2 is a composite of SBS3 and SBS8 (best match SBS8); P3 is SBS3
  pred <- withr::with_seed(9, cbind(
    P1 = noisy_sig(sigs[, "SBS1"]),
    P2 = noisy_sig(rowMeans(sigs[, c("SBS3", "SBS8")]), 0.0005),
    P3 = noisy_sig(sigs[, "SBS3"])
  ))
  pred_exposures <- rbind(
    P1 = truth$exposures["SBS1", ],
    P2 = truth$exposures["SBS3", ] * 0.5 + truth$exposures["SBS8", ],
    P3 = truth$exposures["SBS3", ] * 0.5
  )
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)

  # SBS3 is matched by P2 and P3: P3 is a duplicate, P2 a duplicate and
  # a composite, each counted once
  summary <- method_view_summary(bench)
  expect_identical(summary[["Num Direct Matched"]], 1L)
  expect_identical(summary[["Num Duplicates"]], 1L)
  expect_identical(summary[["Num Dup/Comp"]], 1L)
  expect_identical(summary[["Num Composites"]], 0L)

  # in Final, P3 and the SBS3 piece of P2 are one signature
  final <- benchmark_get_comparison(bench, "m", "Final")
  expect_identical(anyDuplicated(final$y_sig_name), 0L)
  final_sigs <- colnames(signatures(
    final_pred(benchmark_get_entry(bench, "m"))
  ))
  merged <- grep("^Merged Signature", final_sigs, value = TRUE)
  expect_length(merged, 1)
  expect_true(grepl("P3", merged) && grepl("SBS3_like", merged))
  expect_false("P3" %in% final_sigs)

  # ggplot2 warns that the fitted lines leave the plots' fixed axes
  plots <- suppressWarnings(benchmark_plot_composite_exposures(bench, "m"))
  expect_length(plots, 1)
})

test_that("Benchmark functions check full_benchmark and method_id", {
  data(full_benchmark_example)
  data(synthetic_breast_counts)
  musica <- create_musica_from_counts(synthetic_breast_counts, "SBS96")
  expect_error(
    benchmark(list(), musica, "m", "SBS96", make_copy = TRUE),
    regexp = "is not a 'full_benchmark' object"
  )
  expect_error(
    benchmark(full_benchmark_example, list(), "m", "SBS96",
              make_copy = TRUE),
    regexp = "'prediction' must be a 'musica' object"
  )

  # each takes (full_benchmark, method_id, ...)
  calls <- list(
    function(fb, id) benchmark_get_entry(fb, id),
    function(fb, id) benchmark_plot_signatures(fb, id, "Final"),
    function(fb, id) benchmark_get_comparison(fb, id, "Final"),
    function(fb, id) benchmark_plot_comparison(fb, id, "Final"),
    function(fb, id) benchmark_plot_exposures(fb, id, "Final"),
    function(fb, id) benchmark_plot_duplicate_exposures(fb, id),
    function(fb, id) benchmark_plot_composite_exposures(fb, id)
  )
  for (f in calls) {
    expect_error(f(list(), "example_k8"),
                 regexp = "is not a 'full_benchmark' object")
    expect_error(f(full_benchmark_example, "missing_id"),
                 regexp = "not found in")
  }
})

test_that("benchmark(plot = TRUE) prints plots, not lists (#262)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS2")])
  withr::local_pdf(withr::local_tempfile(fileext = ".pdf"))
  plot_output <- function(pred, pred_exposures) {
    setup <- setup_benchmark(truth, pred, pred_exposures)
    utils::capture.output(suppressMessages(suppressWarnings(
      benchmark(setup$bench, setup$prediction, "pred", "SBS96",
        method_id = "m", plot = TRUE, make_copy = TRUE
      )
    )))
  }

  # no duplicates or composites: the plot lists are empty
  exact <- plot_output(truth$signatures, truth$exposures)
  expect_false(any(grepl("list()", exact, fixed = TRUE)))

  # a duplicate: its figure is printed without a list index
  dup <- withr::with_seed(2, noisy_sig(sigs[, "SBS1"]))
  duplicated <- plot_output(
    cbind(truth$signatures, dupA = dup),
    rbind(truth$exposures, dupA = 10)
  )
  expect_false(any(grepl("[[1]]", duplicated, fixed = TRUE)))

  # each duplicate and composite figure is printed once
  local_mocked_bindings(
    benchmark_plot_duplicate_exposures = function(...) list("dup figure"),
    benchmark_plot_composite_exposures = function(...) list("comp figure")
  )
  mocked <- plot_output(truth$signatures, truth$exposures)
  expect_identical(
    grep("figure", mocked, value = TRUE),
    c("[1] \"dup figure\"", "[1] \"comp figure\"")
  )
})

test_that("A duplicate group with no exposure is merged unweighted (#263)", {
  data(cosmic_v2_sigs)
  sigs <- signatures(cosmic_v2_sigs)
  truth <- sim_truth(sigs[, c("SBS1", "SBS2", "SBS3")])
  dup <- withr::with_seed(2, noisy_sig(sigs[, "SBS1"]))
  pred <- cbind(truth$signatures, dupA = dup)
  pred_exposures <- rbind(truth$exposures, dupA = 0)
  pred_exposures["SBS1", ] <- 0
  setup <- setup_benchmark(truth, pred, pred_exposures)
  bench <- run_benchmark(setup$bench, setup$prediction)

  intermediate <- intermediate_pred(benchmark_get_entry(bench, "m"))
  merged <- signatures(intermediate)[, "Merged Signature (SBS1.dupA)"]
  expected <- (sigs[, "SBS1"] + dup) / sum(sigs[, "SBS1"] + dup)
  expect_equal(unname(merged), unname(expected))
  summary <- method_view_summary(bench)
  expect_true(is.finite(summary[["Final RE"]]))
  expect_identical(summary[["Num Missing"]], 0L)
})

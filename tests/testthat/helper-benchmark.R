# Small simulated benchmarks for test-benchmarking.R

# True exposures (signatures x samples) and counts for 'signatures'
sim_truth <- function(signatures, n_samples = 20, seed = 1) {
  samples <- paste0("sample", seq_len(n_samples))
  exposures <- withr::with_seed(seed, matrix(
    sample(200:2000, ncol(signatures) * n_samples, replace = TRUE),
    nrow = ncol(signatures),
    dimnames = list(colnames(signatures), samples)
  ))
  counts <- round(signatures %*% exposures)
  list(signatures = signatures, exposures = exposures, counts = counts)
}

# A signature with small non-negative noise added, renormalised
noisy_sig <- function(signature, sd = 0.002) {
  noisy <- signature + abs(stats::rnorm(length(signature), 0, sd))
  return(noisy / sum(noisy))
}

# A full_benchmark for 'truth' made with 'variant_class', and a musica
# holding the prediction as model 'pred' of 'modality'
setup_benchmark <- function(truth, pred_sigs, pred_exposures,
                            modality = "SBS96", variant_class = modality) {
  bench <- create_benchmark(
    truth$signatures, truth$exposures, truth$counts, variant_class
  )
  prediction <- create_musica_from_counts(truth$counts, variant_class)
  add_result(
    pred_sigs, pred_exposures, prediction, "result", modality, "pred"
  )
  return(list(bench = bench, prediction = prediction))
}

# Benchmark model 'pred' as method "m", without plots or printed output
run_benchmark <- function(bench, prediction, modality = "SBS96", ...) {
  result <- NULL
  utils::capture.output(suppressMessages(
    result <- benchmark(bench, prediction, "pred", modality,
      method_id = "m", plot = FALSE, make_copy = TRUE, ...
    )
  ))
  return(result)
}
